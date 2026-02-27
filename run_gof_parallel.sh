#!/usr/bin/env bash
set -euo pipefail

# =========================
# Config
# =========================
DEFAULT_CARDS_DIR="/afs/cern.ch/user/l/lotan/smqawa/DCTools/cards-eft_ZZ2l2nu"
CARDS_DIR="${CARDS_DIR:-${DEFAULT_CARDS_DIR}}"

OP="${1:-FT8}"
SEED="${2:-123456}"
NTOYS="${3:-500}"

# IMPORTANT with -u: define before cmsset_default
export VO_CMS_SW_DIR="${VO_CMS_SW_DIR:-/cvmfs/cms.cern.ch}"
source /cvmfs/cms.cern.ch/cmsset_default.sh

COMBINE_SETUP="${COMBINE_SETUP:-/afs/cern.ch/user/l/lotan/Combine.sh}"
if [[ -f "${COMBINE_SETUP}" ]]; then
  source "${COMBINE_SETUP}"
else
  echo "[ERROR] Combine setup script not found at ${COMBINE_SETUP}." >&2
  exit 1
fi

# =========================
# Inputs
# =========================
WS="workspace_${OP}_Run2_shape_SR_CR_combined.root"

# 本地 / condor 兼容
if [[ -f "${WS}" ]]; then
  echo ">>> Using local workspace: ${WS}"
else
  cd "${CARDS_DIR}"
  echo ">>> Switched to CARDS_DIR: $(pwd)"
fi

[[ -f "${WS}" ]] || { echo "[ERROR] Missing workspace: ${WS}"; exit 1; }

echo ">>> Working dir: $(pwd)"
echo ">>> Host: $(hostname)"
echo ">>> Start time: $(date)"
echo ">>> combine: $(which combine)"
combine -h | head -n 3 || true

# =========================
# Minimizer / stability options
# =========================
MINOPTS=(
  --cminDefaultMinimizerStrategy=0
  --X-rtd MINIMIZER_analytic
  --X-rtd MINIMIZER_MaxCalls=99999999999
  --cminFallbackAlgo Minuit2,Migrad,0:0.3
  --X-rtd FITTER_NEVER_GIVE_UP
  --X-rtd FITTER_BOUND
)

# =========================
# 1) GoF on data
# =========================
echo "-> GoF (data) ${OP}"
combine -M GoodnessOfFit "${WS}" \
  --algo saturated \
  --redefineSignalPOIs "${OP}" \
  --setParameters r=1 \
  --freezeParameters r \
  -m 125 \
  -n ".gof_data_${OP}" \
  "${MINOPTS[@]}" \
  > "gof_data_${OP}.log" 2>&1

# =========================
# 2) GoF toys
# =========================
echo "-> GoF (toys) ${OP}, ntoys=${NTOYS}, seed=${SEED}"
combine -M GoodnessOfFit "${WS}" \
  --algo saturated \
  --redefineSignalPOIs "${OP}" \
  --setParameters r=1 \
  --freezeParameters r \
  -m 125 \
  -n ".gof_toys_${OP}" \
  -t "${NTOYS}" -s "${SEED}" \
  --toysFrequentist \
  "${MINOPTS[@]}" \
  > "gof_toys_${OP}.log" 2>&1

# =========================
# 3) Collect + plot
# =========================
DATA_ROOT="higgsCombine.gof_data_${OP}.GoodnessOfFit.mH125.root"
TOYS_ROOT="higgsCombine.gof_toys_${OP}.GoodnessOfFit.mH125.${SEED}.root"

[[ -f "${DATA_ROOT}" ]] || { echo "[ERROR] Missing output: ${DATA_ROOT}"; exit 1; }
[[ -f "${TOYS_ROOT}" ]] || { echo "[ERROR] Missing output: ${TOYS_ROOT}"; exit 1; }

echo "-> CollectGoF"
combineTool.py -M CollectGoodnessOfFit \
  --input "${DATA_ROOT}" "${TOYS_ROOT}" \
  -m 125 \
  -o "gof_${OP}.json"

echo "-> plotGof"
plotGof.py "gof_${OP}.json" \
  --statistic saturated \
  --mass 125.0 \
  --title-right="${OP}" \
  --output "gof_${OP}"

echo ">>> Done: gof_${OP}.*"
echo ">>> End time: $(date)"
