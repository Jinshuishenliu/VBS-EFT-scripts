#!/usr/bin/env bash
set -euo pipefail

# Use plot1DScan.py from PATH by default (lxplus/Combine env provides it)
PLOT_TOOL="${PLOT_TOOL:-plot1DScan.py}"

# Default operator list (array!)
DEFAULT_OPS=(FM0 FM1 FM2 FM3 FM4 FM5 FM7 FS0 FS1 FS2 FT0 FT1 FT2 FT3 FT4 FT5 FT6 FT7 FT8 FT9)

# If user passed args -> use them, else use defaults
if [[ $# -gt 0 ]]; then
  OPS=("$@")
else
  OPS=("${DEFAULT_OPS[@]}")
fi

# Ensure plot tool exists
if ! command -v "${PLOT_TOOL}" >/dev/null 2>&1; then
  echo "[ERROR] ${PLOT_TOOL} not found in PATH."
  echo "       Try: export PLOT_TOOL=/full/path/to/plot1DScan.py"
  exit 1
fi

# Resolve actual path (nice for debugging)
PLOT_TOOL_PATH="$(command -v "${PLOT_TOOL}")"
echo ">>> Using plot tool: ${PLOT_TOOL_PATH}"

for OP in "${OPS[@]}"; do
  OUT_FULL="higgsCombine_${OP}_Run2_fullsyst_exp.MultiDimFit.mH125.root"
  OUT_STAT="higgsCombine_${OP}_Run2_statonly_exp.MultiDimFit.mH125.root"
  PLOTOUT="${OP}_shape_syst_vs_stat_exp"

  if [[ ! -f "${OUT_FULL}" ]]; then
    echo "[SKIP] missing ${OUT_FULL}"
    continue
  fi
  if [[ ! -f "${OUT_STAT}" ]]; then
    echo "[SKIP] missing ${OUT_STAT}"
    continue
  fi

  echo ">>> plotting ${OP} -> ${PLOTOUT}.(png|pdf)"

  "${PLOT_TOOL}" "${OUT_FULL}" \
    --POI "${OP}" \
    --main-label "Total exp." \
    --main-color 1 \
    --others "${OUT_STAT}:Stat. only:2" \
    --breakdown syst,stat \
    -o "${PLOTOUT}"
done
