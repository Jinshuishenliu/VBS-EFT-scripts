#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Automatically build Combine datacards for ZZ→2ℓ2ν EFT scans.

The script scans the standard Run-2 years (2016, 2016APV, 2017, 2018), reads the
existing SM datacards to obtain background rates and the observed yield, and
writes new datacards that reference the corresponding AAC templates.
"""

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

OPERATORS: Sequence[str] = (
    "FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7",
    "FS0", "FS1", "FS2",
    "FT0", "FT1", "FT2", "FT3", "FT4", "FT5", "FT6", "FT7", "FT8", "FT9",
)

# OPERATORS: Sequence[str] = (
#     "FT8",
# )

BACKGROUND_PROCS: List[str] = ["WW", "WZ", "DY", "Top", "VVV", "ggVV", "ZZ2l2nu"]

BASE_PROC_IDS = {
    "WW": 1,
    "WZ": 2,
    "DY": 3,
    "Top": 4,
    "VVV": 5,
    "ggVV": 6,
    "ZZ2l2nu": 7,
}

DEFAULT_OBS = 49.0

BASE_DEFAULT_RATES = {
    "WW": 8.541,
    "WZ": 25.538,
    "DY": 62.829,
    "Top": 74.916,
    "VVV": 0.930,
    "ggVV": 4.829,
    "ZZ2l2nu": 17.047,
}

STANDARD_YEARS: Sequence[str] = ("2016", "2016APV", "2017", "2018")
BACKGROUND_NAMES = {"WW", "WZ", "DY", "Top", "VVV", "ggVV", "ZZ2l2nu"}


def build_procs(operator: str) -> List[str]:
    """Return ordered process list for the given operator."""
    return BACKGROUND_PROCS + ["sm", f"sm_lin_quad_{operator}", f"quad_{operator}"]


def build_proc_ids(operator: str) -> Dict[str, int]:
    """Return process ID mapping for Combine (backgrounds positive, signals non-positive)."""
    proc_ids = dict(BASE_PROC_IDS)
    proc_ids.update({
        "sm": 0,
        f"sm_lin_quad_{operator}": -1,
        f"quad_{operator}": -2,
    })
    return proc_ids


def build_default_rates(operator: str) -> Dict[str, float]:
    rates = dict(BASE_DEFAULT_RATES)
    rates.update({
        "sm": -1,
        f"sm_lin_quad_{operator}": -1,
        f"quad_{operator}": -1,
    })
    return rates


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
    signal_procs: Sequence[str],
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
        elif proc in signal_procs:
            applies[proc] = signal_value if signal_value is not None else "-"
        else:
            applies[proc] = "-"
    return applies


def generate_datacard(
    *,
    year: str,
    operator: str,
    bin_name: str,
    output_path: Path,
    aac_root: Path,
    bkg_root: Path,
    rate_source: Optional[Path],
) -> Optional[Path]:
    procs = build_procs(operator)
    proc_ids = build_proc_ids(operator)
    signal_procs = [p for p in procs if p not in BACKGROUND_NAMES]

    rates = build_default_rates(operator)
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
                if proc in rates and proc not in signal_procs:
                    rates[proc] = value
            # if "eft_ZZ2l2nu" in parsed_rates:
            #     rates["sm"] = parsed_rates["eft_ZZ2l2nu"]
        source_processes = processes
        systematics = parsed_systematics
        trailing_lines = rest_lines

    widths = [14] + [16] * len(procs)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "imax * number of categories",
        "jmax * number of samples minus one",
        "kmax * number of nuisance parameters",
        "------------------------------",
        f"# 信号三进程 → aac_{operator}_templates_{year}.root 中的 sm / sm_lin_quad_{operator} / quad_{operator}",
        f"shapes sm                   *  {aac_root}  sm           sm_$SYSTEMATIC",
        f"shapes sm_lin_quad_{operator}          *  {aac_root}  sm_lin_quad_{operator}  sm_lin_quad_{operator}_$SYSTEMATIC",
        f"shapes quad_{operator}                 *  {aac_root}  quad_{operator}         quad_{operator}_$SYSTEMATIC",
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
            applies = map_systematic_values(mapping, source_processes, procs, signal_procs)
            # for sig in ("sm", "sm_lin_quad_FM2", "quad_FM2"):
            #     applies[sig] = "-"
            lines.append(sys_line(name, stype, applies, procs))
    else:
        # Fallback to a minimal placeholder when the source card is missing.
        lines.append(sys_line(f"CMS_SMP23001_Interference_{year}", "lnN", {
            "sm": "1.080", f"sm_lin_quad_{operator}": "1.080", f"quad_{operator}": "1.080",
        }, procs))

    corr_line = f"{bin_name} autoMCCorr aac_{operator}_templates_{year}.root histo_correlation"

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
            f"{bin_name} autoMCStats 0 1 1",
            corr_line,
            "",
        ])

    text = "\n".join(str(line) for line in lines)
    output_path.write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
    print(f"Wrote datacard to {output_path}")
    return output_path


def main() -> None:
    for operator in OPERATORS:
        for year in STANDARD_YEARS:
            bin_name = f"vbs-SR{year}sm"
            output_path = Path(f"shapes-vbs-SR{year}sm_{operator}_shape.dat")
            aac_root = Path(f"aac_{operator}_templates_{year}.root")
            bkg_root = Path(f"shapes-vbs-SR{year}sm.root")
            rate_source = Path(f"shapes-vbs-SR{year}sm.dat")

            if not rate_source.exists():
                print(f"[WARN] {rate_source} 未找到，跳过 {year}/{operator}。")
                continue

            if not aac_root.exists():
                print(f"[WARN] AAC ROOT '{aac_root}' 对应 {year}/{operator} 不存在。")
            if not bkg_root.exists():
                print(f"[WARN] 背景 ROOT '{bkg_root}' 对应 {year}/{operator} 不存在。")

            generate_datacard(
                year=year,
                operator=operator,
                bin_name=bin_name,
                output_path=output_path,
                aac_root=aac_root,
                bkg_root=bkg_root,
                rate_source=rate_source,
            )


if __name__ == "__main__":
    main()
