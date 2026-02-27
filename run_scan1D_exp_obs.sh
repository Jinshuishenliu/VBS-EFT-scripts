#!/usr/bin/env bash
set -euo pipefail

DEFAULT_CARDS_DIR="/afs/cern.ch/user/l/lotan/smqawa/DCTools/cards-eft_ZZ2l2nu"
CARDS_DIR="${CARDS_DIR:-${DEFAULT_CARDS_DIR}}"

OP="${1:-FT8}"
USER_RANGE="${2:-}"     # optional override for RANGE only
STEP="${STEP:-0.1}"

ABS_MAX_DEFAULT="10"
ABS_MAX="${ABS_MAX_DEFAULT}"
case "${OP}" in
  FM0) ABS_MAX="10" ;; FM1) ABS_MAX="30" ;; FM2) ABS_MAX="10" ;; FM3) ABS_MAX="40" ;;
  FM4) ABS_MAX="30" ;; FM5) ABS_MAX="40" ;; FM7) ABS_MAX="60" ;;
  FS0) ABS_MAX="50" ;; FS1) ABS_MAX="40" ;; FS2) ABS_MAX="50" ;;
  FT0) ABS_MAX="2" ;;  FT1) ABS_MAX="2" ;;  FT2) ABS_MAX="5" ;;  FT3) ABS_MAX="5" ;;
  FT4) ABS_MAX="10" ;; FT5) ABS_MAX="5" ;;  FT6) ABS_MAX="8" ;;  FT7) ABS_MAX="15" ;;
  FT8) ABS_MAX="5" ;;  FT9) ABS_MAX="8" ;;
esac

RANGE_DEFAULT="-$ABS_MAX,$ABS_MAX"
RANGE="${USER_RANGE:-${RANGE_DEFAULT}}"

POINTS_OVERRIDE=1000
POINTS="${POINTS_OVERRIDE}"
echo ">>> [INFO] Forcing scan points to ${POINTS}"

export VO_CMS_SW_DIR="${VO_CMS_SW_DIR:-/cvmfs/cms.cern.ch}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
source "${COMBINE_SETUP:-/afs/cern.ch/user/l/lotan/Combine.sh}"

PLOT_TOOL="$(which plot1DScan.py 2>/dev/null || true)"
[[ -n "${PLOT_TOOL}" ]] || { echo "[ERROR] plot1DScan.py not found in PATH"; exit 1; }

MINOPTS=(
  --robustFit 1 --setRobustFitTolerance 0.2
  --cminDefaultMinimizerStrategy 0
  --X-rtd MINIMIZER_analytic
  --X-rtd MINIMIZER_MaxCalls=99999999999
  --cminFallbackAlgo Minuit2,Migrad,0:0.2
  --stepSize 0.005
  --X-rtd FITTER_NEVER_GIVE_UP
  --X-rtd FITTER_BOUND
)

# --- workspace location like your reference script ---
WS="workspace_${OP}_Run2_shape_SR_CR_combined.root"
if [[ -f "${WS}" ]]; then
  echo ">>> Using local workspace in $(pwd)"
else
  cd "${CARDS_DIR}"
  echo ">>> Switched to CARDS_DIR: $(pwd)"
fi
[[ -f "${WS}" ]] || { echo "[ERROR] Missing workspace: ${WS}"; exit 1; }

echo ">>> wd=$(pwd)"
echo ">>> OP=${OP}  range=${RANGE}  points=${POINTS}"
echo ">>> start: $(date)"

# ----------------------------------------------------
# 0) init snapshot for EXPECTED (Asimov)
# ----------------------------------------------------
SNAP="MultiDimFit"
TAG_INIT_EXP="_${OP}_Run2_init_exp"
OUT_INIT_EXP="higgsCombine${TAG_INIT_EXP}.MultiDimFit.mH125.root"

if [[ ! -f "${OUT_INIT_EXP}" ]]; then
  echo ">>> [INIT] expected snapshot: ${SNAP}"
  combine -M MultiDimFit "${WS}" \
    --saveWorkspace \
    --redefineSignalPOIs "${OP}" \
    --setParameters r=1 \
    --freezeParameters r \
    --setParameterRanges "${OP}=${RANGE}" \
    -t -1 -m 125 -n "${TAG_INIT_EXP}" \
    "${MINOPTS[@]}" \
    > "log_init_${OP}_exp.txt" 2>&1
fi
[[ -f "${OUT_INIT_EXP}" ]] || { echo "[ERROR] Missing ${OUT_INIT_EXP}"; exit 1; }

# ----------------------------------------------------
# 1) expected scan from expected snapshot
# ----------------------------------------------------
TAG_SCAN_EXP="_${OP}scan_exp_Run2_SR_CR_combined"
OUT_EXP="higgsCombine${TAG_SCAN_EXP}.MultiDimFit.mH125.root"

echo ">>> [SCAN] expected (Asimov) from snapshot"
combineTool.py "${OUT_INIT_EXP}" -M MultiDimFit \
  --snapshotName "${SNAP}" --skipInitialFit \
  --redefineSignalPOIs "${OP}" \
  --setParameters r=1 \
  --freezeParameters r \
  --algo grid --points "${POINTS}" \
  --setParameterRanges "${OP}=${RANGE}" \
  -t -1 -m 125 -n "${TAG_SCAN_EXP}" \
  "${MINOPTS[@]}" \
  > "log_scan_${OP}_exp.txt" 2>&1

[[ -f "${OUT_EXP}" ]] || { echo "[ERROR] Missing ${OUT_EXP}"; exit 1; }

# ----------------------------------------------------
# 2) init snapshot for OBSERVED (data)
# ----------------------------------------------------
TAG_INIT_OBS="_${OP}_Run2_init_obs"
OUT_INIT_OBS="higgsCombine${TAG_INIT_OBS}.MultiDimFit.mH125.root"

if [[ ! -f "${OUT_INIT_OBS}" ]]; then
  echo ">>> [INIT] observed snapshot: ${SNAP}"
  combine -M MultiDimFit "${WS}" \
    --saveWorkspace \
    --redefineSignalPOIs "${OP}" \
    --setParameters r=1 \
    --freezeParameters r \
    --setParameterRanges "${OP}=${RANGE}" \
    -m 125 -n "${TAG_INIT_OBS}" \
    "${MINOPTS[@]}" \
    > "log_init_${OP}_obs.txt" 2>&1
fi
[[ -f "${OUT_INIT_OBS}" ]] || { echo "[ERROR] Missing ${OUT_INIT_OBS}"; exit 1; }

# ----------------------------------------------------
# 3) observed scan from observed snapshot
# ----------------------------------------------------
TAG_SCAN_OBS="_${OP}scan_obs_Run2_SR_CR_combined"
OUT_OBS="higgsCombine${TAG_SCAN_OBS}.MultiDimFit.mH125.root"

echo ">>> [SCAN] observed (data) from snapshot"
combineTool.py "${OUT_INIT_OBS}" -M MultiDimFit \
  --snapshotName "${SNAP}" --skipInitialFit \
  --redefineSignalPOIs "${OP}" \
  --setParameters r=1 \
  --freezeParameters r \
  --algo grid --points "${POINTS}" \
  --setParameterRanges "${OP}=${RANGE}" \
  -m 125 -n "${TAG_SCAN_OBS}" \
  "${MINOPTS[@]}" \
  > "log_scan_${OP}_obs.txt" 2>&1

[[ -f "${OUT_OBS}" ]] || { echo "[ERROR] Missing ${OUT_OBS}"; exit 1; }

# ----------------------------------------------------
# 4) Plot: overlay exp vs obs
# ----------------------------------------------------
PLOTOUT="${OP}scan_Run2_SR_CR_combined_exp_obs"

python3 "${PLOT_TOOL}" "${OUT_EXP}" \
  --POI "${OP}" \
  --main-label "Expected" \
  --main-color 1 \
  --others "${OUT_OBS}:Observed:2" \
  -o "${PLOTOUT}" \
  > "log_${OP}_exp_vs_obs.txt" 2>&1

echo ">>> Done: ${PLOTOUT}.(png|pdf)"
echo ">>> end: $(date)"
