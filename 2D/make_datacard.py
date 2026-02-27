#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Automatically build Combine datacards for ZZ→2ℓ2ν EFT scans.

The script loops over the standard Run-2 years, looks for AAC template ROOT files
(`aac_<operator token>_templates_<year>.root`), and, for each template, writes a
datacard that references all signal components contained in that template file.

Signal components are built dynamically:
    - sm
    - sm_lin_quad_<operator>, quad_<operator>  (for every operator)
    - mix_<operator_i>_<operator_j>            (for every operator pair)

All background information (rates, systematics) is inherited from the
corresponding SM datacard `shapes-vbs-SR<year>sm.dat`.
"""

from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import re

BACKGROUND_PROCS: List[str] = ["WW", "WZ", "DY", "Top", "VVV", "ggVV", "ZZ2l2nu"]
DEFAULT_OBS = 49.0
DEFAULT_RATES: Dict[str, float] = {
    "WW": 8.541,
    "WZ": 25.538,
    "DY": 62.829,
    "Top": 74.916,
    "VVV": 0.930,
    "ggVV": 4.829,
    "ZZ2l2nu": 17.047,
}
DEFAULT_SIGNAL_RATE = -1.0

STANDARD_YEARS: Sequence[str] = ("2016", "2016APV", "2017", "2018")
BACKGROUND_NAMES = set(BACKGROUND_PROCS)


def fmt_row(cells: Sequence[object], widths: Sequence[int]) -> str:
    return " ".join(str(cell).rjust(width) for cell, width in zip(cells, widths))


def sys_line(name: str, stype: str, applies: Dict[str, str], procs: Sequence[str]) -> str:
    widths = [28, 8] + [8] * len(procs)
    row = [name, stype] + [applies.get(proc, "-") for proc in procs]
    return fmt_row(row, widths)


def extract_datacard_info(
    path: Path,
) -> Tuple[Optional[float], Dict[str, float], Sequence[str], List[Tuple[str, str, Dict[str, str]]], List[str]]:
    """Read observation, rates, systematics, and trailing lines from an SM datacard."""
    observation: Optional[float] = None
    processes: List[str] = []
    rates: Dict[str, float] = {}
    tail_raw: List[str] = []
    rate_section = False

    try:
        with open(path, encoding="utf-8") as src:
            for raw_line in src:
                line = raw_line.rstrip("\n")
                stripped = line.strip()
                if not stripped:
                    if rate_section:
                        tail_raw.append(line)
                    continue

                tokens = stripped.split()

                if observation is None and len(tokens) >= 2 and tokens[0].lower() == "observation":
                    try:
                        observation = float(tokens[1])
                    except ValueError:
                        pass
                    continue

                if not processes and tokens and tokens[0] == "process":
                    if len(tokens) > 1 and not tokens[1].lstrip("-").isdigit():
                        processes = tokens[1:]
                    continue

                if tokens and tokens[0] == "rate":
                    for proc, token in zip(processes, tokens[1:]):
                        try:
                            rates[proc] = float(token)
                        except ValueError:
                            continue
                    rate_section = True
                    continue

                if rate_section:
                    if stripped.startswith("-"):
                        continue
                    tail_raw.append(line)
    except FileNotFoundError:
        return observation, rates, processes, [], []

    systematics: List[Tuple[str, str, Dict[str, str]]] = []
    trailing_lines: List[str] = []
    for line in tail_raw:
        stripped = line.strip()
        if not stripped:
            trailing_lines.append("")
            continue

        tokens = stripped.split()
        if len(tokens) >= 2 and tokens[1] in {"shape", "lnN", "lnU", "gmN", "gmM", "gmG"}:
            mapping = {proc: val for proc, val in zip(processes, tokens[2:])}
            systematics.append((tokens[0], tokens[1], mapping))
        else:
            trailing_lines.append(line)

    return observation, rates, processes, systematics, trailing_lines


def map_systematic_values(
    mapping: Dict[str, str],
    source_processes: Sequence[str],
    procs: Sequence[str],
    signal_names: Sequence[str],
) -> Dict[str, str]:
    """Convert a systematic row from the SM card to the expanded process list."""
    signal_value = None
    for proc in source_processes:
        if proc not in BACKGROUND_NAMES and proc in mapping:
            signal_value = mapping[proc]
            break
    if signal_value is None and "eft_ZZ2l2nu" in mapping:
        signal_value = mapping["eft_ZZ2l2nu"]

    applies: Dict[str, str] = {}
    for proc in procs:
        if proc in mapping:
            applies[proc] = mapping[proc]
        elif proc in signal_names:
            applies[proc] = signal_value if signal_value is not None else "-"
        else:
            applies[proc] = "-"
    return applies


def build_process_lists(operators: Sequence[str]) -> Tuple[List[str], Dict[str, int], List[str]]:
    signal_order: List[str] = ["sm"]
    for op in operators:
        signal_order.append(f"sm_lin_quad_{op}")
        signal_order.append(f"quad_{op}")
    for op_i, op_j in combinations(operators, 2):
        signal_order.append(f"sm_lin_quad_mixed_{op_i}_{op_j}")

    procs = BACKGROUND_PROCS + signal_order
    proc_ids: Dict[str, int] = {}
    for idx, proc in enumerate(BACKGROUND_PROCS, start=1):
        proc_ids[proc] = idx
    proc_ids["sm"] = 0
    next_neg = -1
    for proc in signal_order[1:]:
        proc_ids[proc] = next_neg
        next_neg -= 1
    return procs, proc_ids, signal_order


def parse_operator_token(aac_root: Path) -> Tuple[str, List[str], str]:
    """
    Extract operator token and year from an AAC filename like
    'aac_FT8_FT9_templates_2018.root'.
    """
    m = re.match(r"aac_(.+)_templates_([0-9APV]+)\.root$", aac_root.name)
    if not m:
        raise ValueError(f"Cannot parse operator token/year from {aac_root.name}")
    operator_token = m.group(1)
    year = m.group(2)
    operators = operator_token.split("_")
    return operator_token, operators, year


def generate_datacard(
    *,
    year: str,
    bin_name: str,
    output_path: Path,
    aac_root: Path,
    operator_token: str,
    operators: Sequence[str],
    bkg_root: Path,
    rate_source: Optional[Path],
) -> Optional[Path]:
    procs, proc_ids, signal_order = build_process_lists(operators)
    signal_names = set(signal_order)

    rates: Dict[str, float] = {proc: DEFAULT_RATES.get(proc, DEFAULT_SIGNAL_RATE) for proc in procs}
    for proc in signal_names:
        rates[proc] = DEFAULT_SIGNAL_RATE

    observation = DEFAULT_OBS
    systematics: List[Tuple[str, str, Dict[str, str]]] = []
    trailing_lines: List[str] = []
    source_processes: Sequence[str] = []

    if rate_source:
        parsed_obs, parsed_rates, processes, parsed_systematics, rest_lines = extract_datacard_info(rate_source)
        if parsed_obs is not None:
            observation = parsed_obs
        if parsed_rates:
            for proc, value in parsed_rates.items():
                if proc in rates and proc not in signal_names:
                    rates[proc] = value
            # if "eft_ZZ2l2nu" in parsed_rates:
            #     rates["sm"] = parsed_rates["eft_ZZ2l2nu"]
        source_processes = processes
        systematics = parsed_systematics
        trailing_lines = rest_lines

    widths = [14] + [16] * len(procs)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    comment_signal = f"# 信号进程 → {aac_root.name} 中的 {', '.join(signal_order)}"
    shape_lines = [
        f"shapes {proc.ljust(22)} *  {aac_root}  {proc.ljust(14)} {proc}_$SYSTEMATIC"
        for proc in signal_order
    ]

    lines = [
        "imax * number of categories",
        "jmax * number of samples minus one",
        "kmax * number of nuisance parameters",
        "------------------------------",
        comment_signal,
        *shape_lines,
        "# 其他背景 → 原文件",
        f"shapes *                    *  {bkg_root}  $PROCESS     $PROCESS_$SYSTEMATIC",
        "",
        f"bin                  {bin_name}",
        f"observation          {observation:.1f}",
        "------------------------------",
        "------------------------------",
        fmt_row(["bin"] + [bin_name] * len(procs), widths),
        fmt_row(["process"] + procs, widths),
        fmt_row(["process"] + [proc_ids[p] for p in procs], widths),
        fmt_row(
            ["rate"] + [
                ("{:.3f}".format(rates[p]) if rates[p] >= 0 else "-1")
                for p in procs
            ],
            widths,
        ),
        "------------------------------",
        "# ================= Systematics =================",
    ]

    if systematics and source_processes:
        for name, stype, mapping in systematics:
            applies = map_systematic_values(mapping, source_processes, procs, signal_names)
            lines.append(sys_line(name, stype, applies, procs))
    else:
        applies = {proc: "-" for proc in procs}
        for proc in signal_names:
            applies[proc] = "1.080"
        lines.append(sys_line(f"CMS_SMP23001_Interference_{year}", "lnN", applies, procs))

    corr_line = f"{bin_name} autoMCCorr aac_FT8_FT9_templates_{year}.root histo_correlation"
    
    if trailing_lines:
        if lines[-1] != "":
            lines.append("")
        normalized_trailing: List[str] = []
        inserted_corr = False
        for raw in trailing_lines:
            stripped = raw.strip()
            if stripped:
                tokens = stripped.split()
                if len(tokens) >= 5 and tokens[0] == bin_name and tokens[1] == "autoMCStats":
                    normalized_trailing.append(f"{bin_name} autoMCStats 10 0 1")
                    normalized_trailing.append(corr_line)
                    inserted_corr = True
                    continue
            normalized_trailing.append(raw)
        if not inserted_corr:
            normalized_trailing.append(corr_line)
        lines.extend(normalized_trailing)
    else:
        lines.extend([
            "",
            "# rateParams（给 WZ/WW/Top 背景的自由度）",
            f"CMS_SMP23001_NormWZ_{year}        rateParam  {bin_name}*   WZ     1   [0.1,10]",
            f"CMS_SMP23001_NormWW_{year}        rateParam  {bin_name}*   WW     1   [0.1,10]",
            f"CMS_SMP23001_NormTop_{year}       rateParam  {bin_name}*   Top    1   [0.1,10]",
            "",
            "# per-bin stats",
            f"{bin_name} autoMCStats 0 0 1",
            "",
        ])

    text = "\n".join(str(line) for line in lines)
    output_path.write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
    print(f"[OK] Wrote datacard → {output_path}")
    return output_path


def main() -> None:
    base_dir = Path(".")
    for year in STANDARD_YEARS:
        rate_source = Path(f"shapes-vbs-SR{year}sm.dat")
        bkg_root = Path(f"shapes-vbs-SR{year}sm.root")

        if not rate_source.exists():
            print(f"[WARN] {rate_source} 未找到，跳过 {year}。")
            continue

        aac_files = sorted(base_dir.glob(f"aac_*_templates_{year}.root"))
        if not aac_files:
            print(f"[WARN] 未在 {base_dir} 找到 {year} 的 AAC 模板，跳过。")
            continue

        for aac_root in aac_files:
            operator_token, operators, year_in_file = parse_operator_token(aac_root)
            if year_in_file != year:
                print(f"[WARN] AAC 文件 {aac_root.name} 年份与循环变量不一致（{year_in_file} vs {year}），跳过。")
                continue

            bin_name = f"vbs-SR{year}sm"
            output_name = f"shapes-vbs-SR{year}sm_{operator_token}_shape.dat"
            generate_datacard(
                year=year,
                bin_name=bin_name,
                output_path=Path(output_name),
                aac_root=aac_root,
                operator_token=operator_token,
                operators=operators,
                bkg_root=bkg_root,
                rate_source=rate_source,
            )


if __name__ == "__main__":
    main()
