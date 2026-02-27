#!/usr/bin/env bash
set -euo pipefail

DEFAULT_CARDS_DIR="/afs/cern.ch/user/l/lotan/smqawa/DCTools/cards-eft_ZZ2l2nu_2D"
CARDS_DIR="${CARDS_DIR:-${DEFAULT_CARDS_DIR}}"

# -------------------------
# Fixed 2D POIs
# -------------------------
XVAR="FT8"
YVAR="FT9"

# Optional overrides:
#   arg1: rangeX  (e.g. -15,15)
#   arg2: rangeY  (e.g. -15,15)
#   arg3: points  (e.g. 1681)
USER_RANGE_X="${1:-}"
USER_RANGE_Y="${2:-}"
USER_POINTS="${3:-}"

RANGE_X="${USER_RANGE_X:-"-15,15"}"
RANGE_Y="${USER_RANGE_Y:-"-15,15"}"
POINTS="${USER_POINTS:-3721}"   # 41x41 by default

# -------------------------
# Combine environment
# -------------------------
export VO_CMS_SW_DIR="${VO_CMS_SW_DIR:-/cvmfs/cms.cern.ch}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
source "${COMBINE_SETUP:-/afs/cern.ch/user/l/lotan/Combine.sh}"
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

# -------------------------
# Inputs/Outputs
# -------------------------
MODEL="HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingEFTNegative:analiticAnomalousCouplingEFTNegative"

DC="datacard_Run2_FT8_FT9_shape_SR_CR_combined.dat"
WS="workspace_Run2_FT8_FT9_shape_SR_CR_combined.root"

TAG="_Run2_FT8_FT9_shape_SR_CR_combined_obs"
OUTROOT="higgsCombine${TAG}.MultiDimFit.mH125.root"

# # -------------------------
# # Where to run
# # -------------------------
# # Prefer local if DC exists, else jump to CARDS_DIR
# if [[ -f "${DC}" || -f "${WS}" ]]; then
#   echo ">>> Using local directory: $(pwd)"
# else
#   cd "${CARDS_DIR}"
#   echo ">>> Switched to CARDS_DIR: $(pwd)"
# fi

# [[ -f "${DC}" ]] || { echo "[ERROR] Missing datacard: ${DC}"; exit 1; }

# # -------------------------
# # Grid sanity check
# # -------------------------
# export POINTS
# GRID_N="$(python3 - <<'PY'
# import math, os
# p=int(os.environ["POINTS"])
# n=int(round(math.sqrt(p)))
# print(n)
# PY
# )"
# if [[ "$((GRID_N*GRID_N))" -ne "${POINTS}" ]]; then
#   echo "[WARN] points=${POINTS} is not a perfect square; grid may be awkward. nearest ~ ${GRID_N}x${GRID_N}=$((GRID_N*GRID_N))"
# else
#   echo ">>> [INFO] 2D grid: ${GRID_N} x ${GRID_N} = ${POINTS} points"
# fi

# echo ">>> POIs=${XVAR},${YVAR}"
# echo ">>> ranges: ${XVAR}=${RANGE_X} ; ${YVAR}=${RANGE_Y}"
# echo ">>> start: $(date)"

# # ============================================================
# # 1) Make workspace (only if missing or older than datacard)
# # ============================================================
# make_ws=false
# if [[ ! -f "${WS}" ]]; then
#   make_ws=true
# else
#   # remake if datacard is newer than workspace
#   if [[ "${DC}" -nt "${WS}" ]]; then
#     make_ws=true
#   fi
# fi

# if [[ "${make_ws}" == true ]]; then
#   echo ">>> [WS] building workspace: ${WS}"
#   text2workspace.py "${DC}" \
#     -P "${MODEL}" \
#     --PO "eftOperators=${XVAR},${YVAR}" \
#     --X-pack-asympows \
#     --optimize-simpdf-constraints=cms \
#     --X-optimizeMHDependency=fixed \
#     --X-allow-no-signal \
#     --X-allow-no-background \
#     -o "${WS}" \
#     --verbose 3 \
#     > "log_t2w_${XVAR}_${YVAR}_obs.txt" 2>&1
# else
#   echo ">>> [WS] workspace exists and is up-to-date: ${WS}"
# fi

# [[ -f "${WS}" ]] || { echo "[ERROR] Missing workspace after text2workspace: ${WS}"; exit 1; }

# # ============================================================
# # 2) Run 2D MultiDimFit grid (Asimov, freeze r)
# # ============================================================
# echo ">>> [SCAN] running 2D grid scan -> ${OUTROOT}"
# # combine -M MultiDimFit "${WS}" \
# #   --redefineSignalPOIs "${XVAR},${YVAR}" \
# #   --setParameterRanges "${XVAR}=${RANGE_X}:${YVAR}=${RANGE_Y}" \
# #   -t -1 --setParameters "r=1,${XVAR}=0,${YVAR}=0" --freezeParameters r \
# #   --algo grid --points "${POINTS}" \
# #   -m 125 -n "${TAG}" \
# #   "${MINOPTS[@]}" \
# #   > "log_scan_${XVAR}_${YVAR}${TAG}.txt" 2>&1
# combine -M MultiDimFit "${WS}" \
#   --redefineSignalPOIs "${XVAR},${YVAR}" \
#   --setParameterRanges "${XVAR}=${RANGE_X}:${YVAR}=${RANGE_Y}" \
#   --setParameters "r=1,${XVAR}=0,${YVAR}=0" --freezeParameters r \
#   --algo grid --points "${POINTS}" \
#   -m 125 -n "${TAG}" \
#   "${MINOPTS[@]}" \
#   > "log_scan_${XVAR}_${YVAR}${TAG}_obs.txt" 2>&1

# [[ -f "${OUTROOT}" ]] || { echo "[ERROR] Missing output root: ${OUTROOT}"; exit 1; }

# ============================================================
# 3) Plot 2D scan
# ============================================================
echo ">>> [PLOT] plotting 2D scan"
python3 plot_2D_scan.py \
  --input "${OUTROOT}" \
  --xvar "${XVAR}" \
  --yvar "${YVAR}" \
  --output "${XVAR}_${YVAR}_Run2_obs" \
  --zmax 12 \
  --lumi 137.6
  # > "log_plot_${XVAR}_${YVAR}${TAG}.txt" 2>&1

echo ">>> [OK] workspace: ${WS}"
echo ">>> [OK] scan root:  ${OUTROOT}"
echo ">>> [OK] plot out:   ${XVAR}_${YVAR}_Run2_obs*"
echo ">>> end: $(date)"
