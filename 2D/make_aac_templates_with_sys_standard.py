#!/usr/bin/env python3
import argparse
import glob
import os
import re
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import ROOT as R
import numpy as np

R.gROOT.SetBatch(True)

# ================== 默认配置（可通过 CLI 覆盖） ==================
DEFAULT_TAG = "2018"
DEFAULT_HIST_NAME = "eft_ZZ2l2nu"
DEFAULT_SM_TEMPLATE = "shapes-vbs-SR{tag}sm.root"
DEFAULT_SIG_TEMPLATE = "shapes-vbs-SR{tag}FT*.root"
DEFAULT_OUT_TEMPLATE = "aac_{operator}_templates_{tag}.root"
DEFAULT_FIT_ABS_MAX = 2.0
DEFAULT_PRESET_NAME = "2018_FT8_FT9"
EPSILON = 1e-12

# ✅ 仅修改：为 FT8/FT9 设置不同的默认筛选范围
FIT_ABS_MAX_BY_OP = {
    "FT8": 5.0,
    "FT9": 8.0,
}

# 预先定义的配置，可根据需要补充 / 修改。
CONFIG_PRESETS = {
    "2018_FT8_FT9": {
        "tag": "2018",
        "operators": ["FT8", "FT9"],
        "fit_abs_max": 8.0,
    },
    "2017_FT8_FT9": {
        "tag": "2017",
        "operators": ["FT8", "FT9"],
        "fit_abs_max": 8.0,
    },
    "2016_FT8_FT9": {
        "tag": "2016",
        "operators": ["FT8", "FT9"],
        "fit_abs_max": 8.0,
    },
    "2016APV_FT8_FT9": {
        "tag": "2016APV",
        "operators": ["FT8", "FT9"],
        "fit_abs_max": 8.0,
    },
    "Run2_FT8_FT9": {
        "batch": [
            "2016APV_FT8_FT9",
            "2016_FT8_FT9",
            "2017_FT8_FT9",
            "2018_FT8_FT9",
        ]
    },
}


@dataclass
class SampleInfo:
    coeffs: Dict[str, float]
    source: Optional[str]  # None 表示 SM 文件
    label: str

    def load_arrays(
        self,
        hist_name: str,
        sm_tfile: R.TFile,
        nbins: int,
        nuis: Optional[str] = None,
        ud: Optional[str] = None,
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        target = hist_name if nuis is None else f"{hist_name}_{nuis}{ud}"
        if self.source is None:
            hist = sm_tfile.Get(target)
            if not hist:
                print(f"[WARN] SM histogram {target} not found; skip.")
                return None
            hist = hist.Clone()
            hist.SetDirectory(0)
        else:
            tf = R.TFile.Open(self.source)
            if not tf or tf.IsZombie():
                print(f"[WARN] Cannot open {self.source}; skip.")
                return None
            hist = tf.Get(target)
            if not hist:
                tf.Close()
                return None
            if hist.GetNbinsX() != nbins:
                print(f"[WARN] Bin mismatch in {os.path.basename(self.source)} for {target}; skip.")
                tf.Close()
                return None
            hist = hist.Clone()
            hist.SetDirectory(0)
            tf.Close()

        y, e2 = hist_to_arrays(hist)
        return y, e2


@dataclass
class SampleData:
    coeffs: Dict[str, float]
    y: np.ndarray
    e2: np.ndarray
    label: str


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild AAC templates for VBS EFT samples. Supports one or multiple operators "
            "and produces positive-defined components (SM, SM+Lin+Quad per operator, Quad per operator, "
            "and Mix terms for operator pairs)."
        )
    )
    parser.add_argument(
        "--tag",
        default=None,
        help=f"Era/tag to use (e.g. 2018, 2017, 2016, 2016APV). Default: {DEFAULT_TAG}",
    )
    parser.add_argument(
        "--operators",
        nargs="+",
        help=(
            "Operators to process (e.g. FT9 or FT8 FT9). If omitted, operators are auto-discovered "
            "and handled one-by-one."
        ),
    )
    parser.add_argument(
        "--hist-name",
        default=None,
        help=f"Base histogram name inside ROOT files. Default: {DEFAULT_HIST_NAME}",
    )
    parser.add_argument(
        "--fit-abs-max",
        type=float,
        default=None,
        help=f"Keep samples with |c| <= this for every active operator. Default: {DEFAULT_FIT_ABS_MAX}",
    )
    parser.add_argument(
        "--fit-keep",
        nargs="+",
        type=float,
        default=None,
        help="Extra coupling values (for any operator) to always keep, e.g. --fit-keep -4 4.",
    )
    parser.add_argument(
        "--sm-file",
        default=None,
        help="SM ROOT file template. Accepts {tag}. Default: shapes-vbs-SR{tag}sm.root",
    )
    parser.add_argument(
        "--sig-pattern",
        default=None,
        help=(
            "Signal ROOT glob template. Accepts {tag} and optionally {operator}. "
            "Default: shapes-vbs-SR{tag}FT*.root"
        ),
    )
    parser.add_argument(
        "--out-template",
        default=None,
        help=(
            "Output ROOT filename template with {tag} and {operator}. "
            "Default: aac_{operator}_templates_{tag}.root "
            "(for multiple operators the token becomes e.g. FT8_FT9)."
        ),
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Name of a preset defined in CONFIG_PRESETS at the top of this script.",
    )
    parser.add_argument(
        "--list-configs",
        action="store_true",
        help="List available preset names (CONFIG_PRESETS) and exit.",
    )
    return parser.parse_args()


def hist_to_arrays(h):
    nb = h.GetNbinsX()
    y = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    e2 = np.array([h.GetBinError(i) ** 2 for i in range(1, nb + 1)], dtype=float)
    e2[e2 <= 0] = EPSILON
    return y, e2


def format_sig_pattern(sig_template: str, tag: str, operator_token: str) -> str:
    fmt_args = {"tag": tag}
    if "{operator}" in sig_template:
        fmt_args["operator"] = operator_token
    return sig_template.format(**fmt_args)


def discover_available_operators(tag: str, sig_template: str) -> List[str]:
    pattern = format_sig_pattern(sig_template, tag, operator_token="FT*")
    operators = set()
    for fn in glob.glob(pattern):
        base = os.path.basename(fn)
        ops = re.findall(r"(FT\d+)", base)
        operators.update(ops)
    return sorted(operators)


def parse_value_token(token: str) -> float:
    if not token:
        return 0.0
    sign = -1.0 if token.startswith("m") else 1.0
    if token[0] in "mpMP":
        token = token[1:]
    token = token.replace("p", ".").replace("P", ".")
    if token in {"", "."}:
        return 0.0
    return sign * float(token)


def parse_coeffs_from_filename(
    filename: str,
    operators: Sequence[str],
) -> Optional[Dict[str, float]]:
    base = os.path.basename(filename)
    pairs = re.findall(r"(FT\d+)_([A-Za-z0-9pP]+)", base)
    if not pairs:
        return None
    ops_in_file = {op for op, _ in pairs}
    allowed = set(operators)
    if ops_in_file - allowed:
        return None
    coeffs = {op: 0.0 for op in operators}
    has_nonzero = False
    for op, token in pairs:
        if op in coeffs:
            val = parse_value_token(token)
            coeffs[op] = val
            if abs(val) > 1e-12:
                has_nonzero = True
    return coeffs if has_nonzero else None


def build_sample_infos(
    tag: str,
    operators: Sequence[str],
    sig_template: str,
) -> List[SampleInfo]:
    pattern = format_sig_pattern(sig_template, tag, operator_token="FT*")
    samples: List[SampleInfo] = []
    # SM sample placeholder
    samples.append(
        SampleInfo(
            coeffs={op: 0.0 for op in operators},
            source=None,
            label="SM",
        )
    )
    for fn in sorted(glob.glob(pattern)):
        coeffs = parse_coeffs_from_filename(fn, operators)
        if not coeffs:
            continue
        samples.append(
            SampleInfo(
                coeffs=coeffs,
                source=fn,
                label=os.path.basename(fn),
            )
        )
    return samples


# ✅ 仅修改：加入 per_op_abs_max（按算符独立 abs max）
def sample_pass_filter(
    coeffs: Dict[str, float],
    fit_abs_max: float,
    fit_keep: Optional[Iterable[float]],
    per_op_abs_max: Optional[Dict[str, float]] = None,
) -> bool:
    if fit_abs_max is None and per_op_abs_max is None:
        return True
    keep_vals = list(fit_keep or [])
    for op, val in coeffs.items():
        if abs(val) < 1e-12:
            continue

        # 优先使用按算符独立上限；否则用统一 fit_abs_max
        if per_op_abs_max is not None and op in per_op_abs_max:
            abs_max = float(per_op_abs_max[op])
        else:
            abs_max = float(fit_abs_max) if fit_abs_max is not None else None

        if abs_max is None:
            continue

        if abs(val) <= abs_max:
            continue

        if any(abs(val - k) < 1e-9 for k in keep_vals):
            continue

        return False
    return True


def gather_sample_data(
    sample_infos: Sequence[SampleInfo],
    hist_name: str,
    sm_tfile: R.TFile,
    nbins: int,
    fit_abs_max: float,
    fit_keep: Optional[Iterable[float]],
    sm_arrays: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    nuis: Optional[str] = None,
    ud: Optional[str] = None,
    per_op_abs_max: Optional[Dict[str, float]] = None,  # ✅ 新增
) -> List[SampleData]:
    data: List[SampleData] = []
    for info in sample_infos:
        if info.source is None and nuis is None and sm_arrays is not None:
            y, e2 = sm_arrays
        else:
            loaded = info.load_arrays(hist_name, sm_tfile, nbins, nuis=nuis, ud=ud)
            if loaded is None:
                continue
            y, e2 = loaded
        if info.source is not None:
            if not sample_pass_filter(info.coeffs, fit_abs_max, fit_keep, per_op_abs_max=per_op_abs_max):
                continue
        data.append(
            SampleData(
                coeffs=dict(info.coeffs),
                y=y.copy(),
                e2=e2.copy(),
                label=info.label,
            )
        )
    return data


def build_design_matrix(
    data: Sequence[SampleData],
    operators: Sequence[str],
    pair_ops: Sequence[Tuple[str, str]],
) -> np.ndarray:
    rows = []
    for dp in data:
        row = [1.0]
        for op in operators:
            row.append(dp.coeffs[op])
        for op in operators:
            row.append(dp.coeffs[op] ** 2)
        for op_i, op_j in pair_ops:
            row.append(dp.coeffs[op_i] * dp.coeffs[op_j])
        rows.append(row)
    return np.array(rows, dtype=float)


def fit_coefficients_per_bin(
    data: Sequence[SampleData],
    operators: Sequence[str],
    pair_ops: Sequence[Tuple[str, str]],
) -> Tuple[
    Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[Tuple[str, str], np.ndarray]],
    List[np.ndarray],
]:
    if not data:
        raise RuntimeError("No sample data provided for fit.")
    nbins = data[0].y.shape[0]
    for dp in data:
        if dp.y.shape[0] != nbins:
            raise RuntimeError("Histogram bin mismatch across samples.")

    design = build_design_matrix(data, operators, pair_ops)
    n_features = design.shape[1]
    if design.shape[0] < n_features:
        raise RuntimeError(
            f"Insufficient sample points ({design.shape[0]}) for {n_features} fit parameters."
        )

    a = np.zeros(nbins)
    lin = {op: np.zeros(nbins) for op in operators}
    quad = {op: np.zeros(nbins) for op in operators}
    mix = {pair: np.zeros(nbins) for pair in pair_ops}
    cov_mats: List[np.ndarray] = []

    for ib in range(nbins):
        y_vec = np.array([dp.y[ib] for dp in data], dtype=float)
        e2_vec = np.array([dp.e2[ib] for dp in data], dtype=float)
        w = 1.0 / e2_vec
        AtWA = design.T @ (w[:, None] * design)
        AtWy = design.T @ (w * y_vec)
        try:
            coeff = np.linalg.solve(AtWA, AtWy)
        except np.linalg.LinAlgError:
            coeff = np.linalg.lstsq(AtWA, AtWy, rcond=None)[0]
        try:
            cov = np.linalg.inv(AtWA)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(AtWA)

        idx = 0
        a[ib] = coeff[idx]
        idx += 1
        for op in operators:
            lin[op][ib] = coeff[idx]
            idx += 1
        for op in operators:
            quad[op][ib] = coeff[idx]
            idx += 1
        for pair in pair_ops:
            mix[pair][ib] = coeff[idx]
            idx += 1
        cov_mats.append(cov)

    # 数值保护：极小负的二次项置零
    for op in operators:
        mask = (quad[op] < 0) & (quad[op] > -1e-10)
        quad[op][mask] = 0.0
    return (a, lin, quad, mix), cov_mats


def write_hist(proto_hist, name: str, values: np.ndarray, tfile_out: R.TFile):
    hist = proto_hist.Clone(name)
    hist.SetDirectory(0)
    hist.Reset()
    for ib, val in enumerate(values, start=1):
        hist.SetBinContent(ib, float(val))
        hist.SetBinError(ib, 0.0)
    tfile_out.cd()
    hist.Write()


def write_positive_components(
    proto_hist,
    operators: Sequence[str],
    pair_ops: Sequence[Tuple[str, str]],
    coeffs: Tuple[np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[Tuple[str, str], np.ndarray]],
    tfile_out: R.TFile,
    suffix: str = "",
):
    a, lin, quad, mix = coeffs
    write_hist(proto_hist, f"sm{suffix}", a, tfile_out)

    single_op = len(operators) == 1
    for op in operators:
        quad_vals = quad[op]
        name_quad = "quad" if single_op else f"quad_{op}"
        write_hist(proto_hist, f"{name_quad}{suffix}", quad_vals, tfile_out)

        slq_vals = a + lin[op] + quad[op]
        name_slq = "sm_lin_quad" if single_op else f"sm_lin_quad_{op}"
        write_hist(proto_hist, f"{name_slq}{suffix}", slq_vals, tfile_out)

    for op_i, op_j in pair_ops:
        mix_vals = a + lin[op_i] + quad[op_i] + lin[op_j] + quad[op_j] + mix[(op_i, op_j)]
        name_mix = f"sm_lin_quad_mixed_{op_i}_{op_j}"
        write_hist(proto_hist, f"{name_mix}{suffix}", mix_vals, tfile_out)


def build_parameter_indices(
    operators: Sequence[str],
    pair_ops: Sequence[Tuple[str, str]],
) -> Tuple[Dict[str, object], int]:
    lin_idx: Dict[str, int] = {}
    quad_idx: Dict[str, int] = {}
    mix_idx: Dict[Tuple[str, str], int] = {}
    current = 0
    sm_idx = current
    current += 1
    for op in operators:
        lin_idx[op] = current
        current += 1
    for op in operators:
        quad_idx[op] = current
        current += 1
    for pair in pair_ops:
        mix_idx[pair] = current
        current += 1
    index_maps: Dict[str, object] = {
        "sm": sm_idx,
        "lin": lin_idx,
        "quad": quad_idx,
        "mix": mix_idx,
    }
    return index_maps, current


def build_template_vectors(
    operators: Sequence[str],
    pair_ops: Sequence[Tuple[str, str]],
) -> Tuple[List[Tuple[str, np.ndarray]], int]:
    index_maps, size = build_parameter_indices(operators, pair_ops)
    single_op = len(operators) == 1

    def make_vec(indices: Iterable[int]) -> np.ndarray:
        vec = np.zeros(size, dtype=float)
        for idx in indices:
            vec[idx] += 1.0
        return vec

    templates: List[Tuple[str, np.ndarray]] = []
    templates.append(("sm", make_vec([index_maps["sm"]])))

    for op in operators:
        quad_name = "quad" if single_op else f"quad_{op}"
        templates.append((quad_name, make_vec([index_maps["quad"][op]])))

        slq_name = "sm_lin_quad" if single_op else f"sm_lin_quad_{op}"
        templates.append(
            (
                slq_name,
                make_vec(
                    [
                        index_maps["sm"],
                        index_maps["lin"][op],
                        index_maps["quad"][op],
                    ]
                ),
            )
        )

    for pair in pair_ops:
        op_i, op_j = pair
        mix_name = f"sm_lin_quad_mixed_{op_i}_{op_j}"
        templates.append(
            (
                mix_name,
                make_vec(
                    [
                        index_maps["sm"],
                        index_maps["lin"][op_i],
                        index_maps["quad"][op_i],
                        index_maps["lin"][op_j],
                        index_maps["quad"][op_j],
                        index_maps["mix"][pair],
                    ]
                ),
            )
        )

    return templates, size


def write_correlation_hist(
    proto_hist,
    operators: Sequence[str],
    pair_ops: Sequence[Tuple[str, str]],
    cov_mats: Sequence[np.ndarray],
    tfile_out: R.TFile,
) -> None:
    templates, size = build_template_vectors(operators, pair_ops)
    nbins = proto_hist.GetNbinsX()
    ntemps = len(templates)
    if ntemps == 0:
        return

    h_corr = R.TH3D(
        "histo_correlation",
        "correlated mc stat uncertainty; xbin; procY; procZ",
        nbins,
        0.5,
        nbins + 0.5,
        ntemps,
        0.5,
        ntemps + 0.5,
        ntemps,
        0.5,
        ntemps + 0.5,
    )

    for idx, (name, _) in enumerate(templates, start=1):
        h_corr.GetYaxis().SetBinLabel(idx, name)
        h_corr.GetZaxis().SetBinLabel(idx, name)

    for ib in range(1, nbins + 1):
        cov = cov_mats[ib - 1]
        if cov.shape[0] != size:
            raise RuntimeError("Covariance matrix size mismatch when building correlation histogram.")
        for iy, (_, vec_y) in enumerate(templates, start=1):
            for iz in range(iy, ntemps + 1):
                vec_z = templates[iz - 1][1]
                val = float(vec_y @ cov @ vec_z)
                h_corr.SetBinContent(ib, iy, iz, val)
                if iz != iy:
                    h_corr.SetBinContent(ib, iz, iy, val)

    tfile_out.cd()
    h_corr.Write()


def describe_samples(data: Sequence[SampleData], operators: Sequence[str]) -> List[str]:
    desc = []
    for dp in data:
        coeff_parts = [f"{op}={dp.coeffs[op]:.6g}" for op in operators]
        desc.append(f"{dp.label}: " + ", ".join(coeff_parts))
    return desc


def run_operator_set(
    tag: str,
    operators: Sequence[str],
    hist_name: str,
    fit_abs_max: float,
    fit_keep: Optional[Iterable[float]],
    sm_template: str,
    sig_template: str,
    out_template: str,
) -> bool:
    operators = sorted(dict.fromkeys(operators))
    if not operators:
        print("[WARN] No operators specified; skip.")
        return False

    sm_file = sm_template.format(tag=tag, operator=operators[0] if operators else "")
    if not os.path.exists(sm_file):
        raise RuntimeError(f"SM file not found: {sm_file}")

    sample_infos = build_sample_infos(tag, operators, sig_template)
    if len(sample_infos) <= 1:
        print(f"[WARN] No signal ROOT files found for operators {operators}; skip.")
        return False

    f_sm = R.TFile.Open(sm_file)
    if not f_sm or f_sm.IsZombie():
        raise RuntimeError(f"Cannot open {sm_file}")

    h_sm_nom = f_sm.Get(hist_name)
    if not h_sm_nom:
        f_sm.Close()
        raise RuntimeError(f"Histogram {hist_name} not found in {sm_file}")

    y_sm_nom, e2_sm_nom = hist_to_arrays(h_sm_nom)
    nbins = h_sm_nom.GetNbinsX()

    # ✅ 仅修改：按算符独立 abs max（FT8=5, FT9=8），其他算符 fallback 到 fit_abs_max
    per_op_abs_max = {op: FIT_ABS_MAX_BY_OP.get(op, fit_abs_max) for op in operators}

    data_nominal = gather_sample_data(
        sample_infos,
        hist_name,
        f_sm,
        nbins,
        fit_abs_max=fit_abs_max,
        fit_keep=fit_keep,
        sm_arrays=(y_sm_nom, e2_sm_nom),
        per_op_abs_max=per_op_abs_max,
    )
    if len(data_nominal) <= 1:
        f_sm.Close()
        raise RuntimeError("Not enough samples after selection to perform the fit.")

    pair_ops = list(combinations(operators, 2))
    # Check presence of mix-sensitive samples if needed
    if pair_ops:
        has_mix_point = any(
            sum(abs(dp.coeffs[op]) > 1e-12 for op in operators) >= 2 for dp in data_nominal
        )
        if not has_mix_point:
            f_sm.Close()
            raise RuntimeError(
                "No samples with multiple operators non-zero were found; cannot extract mix terms."
            )

    print("[INFO] Fit (nominal) uses coefficients:")
    for line in describe_samples(data_nominal, operators):
        print(f"       - {line}")

    coeffs_nom, cov_mats = fit_coefficients_per_bin(data_nominal, operators, pair_ops)

    operator_token = "_".join(operators)
    out_file = out_template.format(tag=tag, operator=operator_token)
    fout = R.TFile(out_file, "RECREATE")

    write_positive_components(
        h_sm_nom,
        operators,
        pair_ops,
        coeffs_nom,
        fout,
        suffix="",
    )
    write_correlation_hist(
        h_sm_nom,
        operators,
        pair_ops,
        cov_mats,
        fout,
    )
    print("[OK] Wrote nominal components:", end=" ")
    parts = ["sm"]
    if len(operators) == 1:
        parts += ["quad", "sm_lin_quad"]
    else:
        for op in operators:
            parts += [f"quad_{op}", f"sm_lin_quad_{op}"]
        for op_i, op_j in pair_ops:
            parts.append(f"sm_lin_quad_mixed_{op_i}_{op_j}")
    print(", ".join(parts))

    nuis_list = list_nuisances_from_sm(f_sm, hist_name)
    print(f"[INFO] Found {len(nuis_list)} nuisances in SM file.")

    for nuis in nuis_list:
        for ud in ("Up", "Down"):
            data_var = gather_sample_data(
                sample_infos,
                hist_name,
                f_sm,
                nbins,
                fit_abs_max=fit_abs_max,
                fit_keep=fit_keep,
                sm_arrays=None,
                nuis=nuis,
                ud=ud,
                per_op_abs_max=per_op_abs_max,  # ✅ 同样应用到系统误差样本筛选
            )
            if len(data_var) <= 1:
                print(f"[WARN] NUIS={nuis}{ud}: not enough samples with this variation; skip.")
                continue
            try:
                coeffs_var, _ = fit_coefficients_per_bin(data_var, operators, pair_ops)
            except RuntimeError as err:
                print(f"[WARN] NUIS={nuis}{ud}: {err}; skip.")
                continue

            suffix = f"_{nuis}{ud}"
            write_positive_components(
                h_sm_nom,
                operators,
                pair_ops,
                coeffs_var,
                fout,
                suffix=suffix,
            )
            print(
                f"[OK] Wrote {nuis}{ud} variations for: "
                + (", ".join(parts))
            )

    fout.Close()
    f_sm.Close()
    print(f"[DONE] Output => {out_file}")
    return True


def list_nuisances_from_sm(sm_tfile, hist_base):
    nuis = set()
    for key in sm_tfile.GetListOfKeys():
        name = key.GetName()
        m = re.match(rf"^{re.escape(hist_base)}_(.+)(Up|Down)$", name)
        if m:
            nuis.add(m.group(1))
    return sorted(nuis)


def main():
    args = parse_args()
    if args.list_configs:
        if not CONFIG_PRESETS:
            print("No presets are defined.")
        else:
            print("Available presets:")
            for name, cfg in CONFIG_PRESETS.items():
                pretty = ", ".join(f"{k}={v}" for k, v in cfg.items())
                print(f"  - {name}: {pretty}")
        return

    # Determine which preset(s) to run
    preset_name = None
    preset = {}
    if args.config:
        if args.config not in CONFIG_PRESETS:
            raise RuntimeError(
                f"Unknown config preset '{args.config}'. Use --list-configs to see available names."
            )
        preset_name = args.config
        preset = CONFIG_PRESETS[preset_name]
    elif DEFAULT_PRESET_NAME and DEFAULT_PRESET_NAME in CONFIG_PRESETS:
        preset_name = DEFAULT_PRESET_NAME
        preset = CONFIG_PRESETS[DEFAULT_PRESET_NAME]

    # Build list of job configurations (each job corresponds to one tag/operator set)
    job_presets: List[Dict[str, object]] = []
    if preset and "batch" in preset:
        for name in preset["batch"]:
            if isinstance(name, str):
                sub = CONFIG_PRESETS.get(name)
                if not sub:
                    raise RuntimeError(f"Preset '{name}' referenced by '{preset_name}' is not defined.")
                job_presets.append(dict(sub))
            elif isinstance(name, dict):
                job_presets.append(dict(name))
            else:
                raise RuntimeError(f"Unsupported entry in batch preset '{preset_name}': {name}")
    else:
        job_presets.append(dict(preset))

    if not job_presets:
        job_presets.append({})

    # Gather CLI overrides (applied to every job)
    override_tag = args.tag
    override_hist = args.hist_name
    override_sig = args.sig_pattern
    override_sm = args.sm_file
    override_out = args.out_template
    override_fit_abs_max = args.fit_abs_max
    override_fit_keep = args.fit_keep
    override_ops = list(args.operators) if args.operators else None

    processed_any = False
    for job in job_presets:
        tag = override_tag if override_tag is not None else job.get("tag", DEFAULT_TAG)
        hist_name = override_hist if override_hist is not None else job.get("hist_name", DEFAULT_HIST_NAME)
        sig_template = (
            override_sig if override_sig is not None else job.get("sig_template", DEFAULT_SIG_TEMPLATE)
        )
        sm_template = override_sm if override_sm is not None else job.get("sm_template", DEFAULT_SM_TEMPLATE)
        out_template = (
            override_out if override_out is not None else job.get("out_template", DEFAULT_OUT_TEMPLATE)
        )

        fit_abs_max = (
            override_fit_abs_max if override_fit_abs_max is not None else job.get("fit_abs_max", DEFAULT_FIT_ABS_MAX)
        )
        fit_keep_values = override_fit_keep if override_fit_keep is not None else job.get("fit_keep")
        fit_keep = set(fit_keep_values or [])

        if override_ops:
            operator_sets = [list(override_ops)]
        else:
            job_ops = job.get("operators")
            if job_ops:
                operator_sets = [list(job_ops)]
            else:
                detected = discover_available_operators(tag, sig_template)
                if not detected:
                    print(
                        f"[WARN] No operators found for tag {tag} using pattern {sig_template}; skip."
                    )
                    continue
                operator_sets = [[op] for op in detected]
                print("[INFO] Auto-discovered operators: " + ", ".join(detected))

        for ops in operator_sets:
            print(f"[INFO] Processing tag={tag} operators={ops}")
            try:
                ok = run_operator_set(
                    tag=tag,
                    operators=ops,
                    hist_name=hist_name,
                    fit_abs_max=fit_abs_max,
                    fit_keep=fit_keep,
                    sm_template=sm_template,
                    sig_template=sig_template,
                    out_template=out_template,
                )
            except RuntimeError as err:
                print(f"[WARN] {tag} operators={ops}: {err}")
                ok = False
            processed_any = processed_any or ok

    if not processed_any:
        raise RuntimeError("No operators were processed successfully.")


if __name__ == "__main__":
    main()

