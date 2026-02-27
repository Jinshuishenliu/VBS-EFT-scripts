#!/usr/bin/env bash
set -euo pipefail

# =========================
# Config
# =========================
DEFAULT_CARDS_DIR="/afs/cern.ch/user/l/lotan/smqawa/DCTools/cards-eft_ZZ2l2nu"
CARDS_DIR="${CARDS_DIR:-${DEFAULT_CARDS_DIR}}"

OP="${1:-FT8}"
PARALLEL="${2:-8}"
USER_RANGE="${3:-}"

# Default POI ranges
RANGE_DEFAULT="-3,3"
case "${OP}" in
  FM0|FM1|FM2|FM3|FM4) RANGE_DEFAULT="-30,30" ;;
  FM5)                 RANGE_DEFAULT="-50,50" ;;
  FM7)                 RANGE_DEFAULT="-100,100" ;;
  FS0|FS1|FS2)          RANGE_DEFAULT="-100,100" ;;
  FT0|FT1)             RANGE_DEFAULT="-8,8" ;;
  FT2|FT3|FT5|FT6|FT8) RANGE_DEFAULT="-10,10" ;;
  FT4|FT9)             RANGE_DEFAULT="-15,15" ;;
  FT7)                 RANGE_DEFAULT="-20,20" ;;
esac
RANGE="${USER_RANGE:-${RANGE_DEFAULT}}"

# Env (keep, since combine needs it)
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

# Condor transfer mode: WS in scratch (CWD), don't cd
# Local mode: cd to CARDS_DIR if WS isn't in CWD
if [[ -f "${WS}" ]]; then
  echo ">>> Using local workspace in $(pwd)"
else
  cd "${CARDS_DIR}"
  echo ">>> Switched to CARDS_DIR: $(pwd)"
fi

echo ">>> Working dir: $(pwd)"
echo ">>> Host: $(hostname)"
echo ">>> Start time: $(date)"
echo ">>> combine: $(which combine)"
combine -h | head -n 3 || true

[[ -f "${WS}" ]] || { echo "[ERROR] Missing workspace: ${WS}"; exit 1; }
echo ">>> Workspace OK: ${WS}"

# =========================
# Impacts (Asimov)
# =========================
MINOPTS="\
  --robustFit 1 --setRobustFitTolerance 0.2 \
  --cminDefaultMinimizerStrategy 0 \
  --X-rtd MINIMIZER_analytic \
  --X-rtd MINIMIZER_MaxCalls=99999999999 \
  --cminFallbackAlgo Minuit2,Migrad,0:0.2 \
  --stepSize 0.005 \
  --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND \
"

# COMMON="\
#   -d ${WS} -m 125 \
#   --redefineSignalPOIs ${OP} \
#   --setParameterRanges ${OP}=${RANGE} \
#   --setParameters r=1,${OP}=0 \
#   --freezeParameters r \
#   -t -1\
#   ${MINOPTS} \
# "
COMMON="\
  -d ${WS} -m 125 \
  --redefineSignalPOIs ${OP} \
  --setParameterRanges ${OP}=${RANGE} \
  --setParameters r=1,${OP}=0 \
  --freezeParameters r \
  ${MINOPTS} \
"

TAG=".impacts_${OP}"

echo ">>> [1/4] Initial fit (Asimov) for impacts"
combineTool.py -M Impacts ${COMMON} --exclude r --doInitialFit -n ${TAG} | tee "log_initialFit_${OP}.txt"

echo ">>> [2/4] Per-nuisance fits (parallel=${PARALLEL})"
combineTool.py -M Impacts ${COMMON} --exclude r --doFits --parallel "${PARALLEL}" -n ${TAG} | tee "log_doFits_${OP}.txt"

echo ">>> [3/4] Collect impacts json"
combineTool.py -M Impacts ${COMMON} --exclude r -o "impacts_${OP}_obs.json" -n ${TAG}

echo ">>> [4/4] Plot"
plotImpacts.py -i "impacts_${OP}_obs.json" -o "impacts_${OP}_obs" --POI "${OP}" --blind

echo ">>> Done: impacts_${OP}.*"
echo ">>> End time: $(date)"
