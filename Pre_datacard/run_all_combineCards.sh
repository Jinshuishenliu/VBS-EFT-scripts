#!/usr/bin/env bash

set -euo pipefail

OPERATORS=(
  "FM0" "FM1" "FM2" "FM3" "FM4" "FM5" "FM7"
  "FS0" "FS1" "FS2"
  "FT0" "FT1" "FT2" "FT3" "FT4" "FT5" "FT6" "FT7"
  "FT8" "FT9"
)

ERAS=("2018" "2017" "2016" "2016APV")

for operator in "${OPERATORS[@]}"; do
  echo "========== ${operator} =========="
  for era in "${ERAS[@]}"; do
    sr_file="shapes-vbs-SR${era}sm_${operator}_shape.dat"
    three_l_file="shapes-vbs-3L${era}.dat"
    em_file="shapes-vbs-EM${era}.dat"
    output_file="shapes-vbs-SR_${era}sm_${operator}_shape_SR_CR_combined.dat"

    missing_files=0
    for f in "$sr_file" "$three_l_file" "$em_file"; do
      if [[ ! -f "$f" ]]; then
        echo "[WARN] Missing input ${f}; skipping ${operator} (${era})" >&2
        missing_files=1
      fi
    done
    (( missing_files )) && continue

    echo "Combining ${operator} for ${era} -> ${output_file}"
    combineCards.py \
      "vbs-SR${era}=${sr_file}" \
      "vbs-3L${era}=${three_l_file}" \
      "vbs-EM${era}=${em_file}" \
      &> "${output_file}"
  done

  run2_output="datacard_Run2_${operator}_shape_SR_CR_combined.dat"
  missing_combined=0
  combine_args=()
  for era in "${ERAS[@]}"; do
    sr_cr_file="shapes-vbs-SR_${era}sm_${operator}_shape_SR_CR_combined.dat"
    if [[ ! -f "$sr_cr_file" ]]; then
      echo "[WARN] Missing combined SR/CR card ${sr_cr_file}; skipping Run2 card for ${operator}" >&2
      missing_combined=1
    fi
    combine_args+=("ch${era}=${sr_cr_file}")
  done

  if (( ! missing_combined )); then
    echo "Combining Run2 card for ${operator} -> ${run2_output}"
    combineCards.py "${combine_args[@]}" > "${run2_output}"
  fi
done

echo "All requested combinations finished."
