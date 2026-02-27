#!/usr/bin/env bash
set -euo pipefail

OPS=(
  FM0 FM1 FM2 FM3 FM4 FM5 FM7
  FS0 FS1 FS2
  FT0 FT1 FT2 FT3 FT4 FT5 FT6 FT7 FT8 FT9
)

MODEL="HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingEFTNegative:analiticAnomalousCouplingEFTNegative"

for OP in "${OPS[@]}"; do
  DC="datacard_Run2_${OP}_shape_SR_CR_combined.dat"
  WS="workspace_${OP}_Run2_shape_SR_CR_combined.root"

  [[ -f "${DC}" ]] || { echo "[WARN] missing ${DC}, skip"; continue; }

  echo "=== ${OP} ==="
  text2workspace.py "${DC}" \
    -P "${MODEL}" \
    --PO "eftOperators=${OP}" \
    --X-pack-asympows \
    --optimize-simpdf-constraints=cms \
    --X-optimizeMHDependency=fixed \
    --X-allow-no-signal \
    --X-allow-no-background \
    -o "${WS}" \
    --verbose 3
done
