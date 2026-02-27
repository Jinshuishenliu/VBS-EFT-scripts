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

# ---- compute POINTS from RANGE and STEP (force to 1000) ----
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

WS="workspace_${OP}_Run2_shape_SR_CR_combined.root"
if [[ -f "${WS}" ]]; then
  echo ">>> Using local workspace in $(pwd)"
else
  cd "${CARDS_DIR}"
  echo ">>> Switched to CARDS_DIR: $(pwd)"
fi
[[ -f "${WS}" ]] || { echo "[ERROR] Missing workspace: ${WS}"; exit 1; }

echo ">>> wd=$(pwd)"
echo ">>> OP=${OP}  range=${RANGE}  step=${STEP}  points=${POINTS}  (EXPECTED Asimov)"
echo ">>> start: $(date)"

# ----------------------------------------------------
# 0) One shared initial fit -> produces ONE snapshot
# ----------------------------------------------------
SNAP="MultiDimFit"
TAG_INIT="_${OP}_Run2_init_exp"
OUT_INIT="higgsCombine${TAG_INIT}.MultiDimFit.mH125.root"

if [[ ! -f "${OUT_INIT}" ]]; then
  echo ">>> [INIT] creating shared snapshot: ${SNAP}"
  combine -M MultiDimFit "${WS}" \
    --saveWorkspace \
    --redefineSignalPOIs "${OP}" \
    --setParameters r=1 \
    --freezeParameters r \
    --setParameterRanges "${OP}=${RANGE}" \
    -t -1 -m 125 -n "${TAG_INIT}" \
    "${MINOPTS[@]}" \
    > "log_init_${OP}_exp.txt" 2>&1
fi
[[ -f "${OUT_INIT}" ]] || { echo "[ERROR] Missing ${OUT_INIT}"; exit 1; }

# helper: run scan from the SAME snapshot
run_scan_from_snap () {
  local tag="$1"          # e.g. _FT8_Run2_fullsyst_exp
  local freeze_extra="$2" # e.g. "" or ",allConstrainedNuisances"
  local logfile="$3"

  combineTool.py "${OUT_INIT}" -M MultiDimFit \
    --snapshotName "${SNAP}" --skipInitialFit \
    --redefineSignalPOIs "${OP}" \
    --setParameters r=1 \
    --freezeParameters "r${freeze_extra}" \
    --algo grid --points "${POINTS}" \
    --setParameterRanges "${OP}=${RANGE}" \
    -t -1 -m 125 -n "${tag}" \
    "${MINOPTS[@]}" \
    > "${logfile}" 2>&1
}

# -------------------------
# 1) FULL syst (expected)  from SAME snapshot
# -------------------------
TAG_FULL="_${OP}_Run2_fullsyst_exp"
OUT_FULL="higgsCombine${TAG_FULL}.MultiDimFit.mH125.root"
echo ">>> [SCAN] full (total) from snapshot"
run_scan_from_snap "${TAG_FULL}" "" "log_scan_${OP}_fullsyst_exp.txt"
[[ -f "${OUT_FULL}" ]] || { echo "[ERROR] Missing ${OUT_FULL}"; exit 1; }

# -------------------------
# 2) STAT only (expected)  from SAME snapshot
# -------------------------
TAG_STAT="_${OP}_Run2_statonly_exp"
OUT_STAT="higgsCombine${TAG_STAT}.MultiDimFit.mH125.root"
echo ">>> [SCAN] stat-only from snapshot"
run_scan_from_snap "${TAG_STAT}" ",allConstrainedNuisances" "log_scan_${OP}_statonly_exp.txt"
[[ -f "${OUT_STAT}" ]] || { echo "[ERROR] Missing ${OUT_STAT}"; exit 1; }

# -------------------------
# 3) Plot: overlay + breakdown
# -------------------------
PLOTOUT="${OP}_shape_syst_vs_stat_exp"

python3 "${PLOT_TOOL}" "${OUT_FULL}" \
  --POI "${OP}" \
  --main-label "Total exp." \
  --main-color 1 \
  --others "${OUT_STAT}:Stat. only:2" \
  --breakdown syst,stat \
  -o "${PLOTOUT}" \
  > "log_plot_${OP}_syst_vs_stat_exp.txt" 2>&1

echo ">>> Done: ${PLOTOUT}.(png|pdf)"
echo ">>> end: $(date)"
