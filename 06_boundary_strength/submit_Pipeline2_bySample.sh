#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/research/plate_c/FINAL"
INPUT_DIR="${PROJECT}/data/merged/input"
SCRIPT="Pipeline2_bySample.sh"

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: cannot find $SCRIPT in $(pwd). cd to the folder containing it, or edit SCRIPT path in this file." >&2
  exit 1
fi

N=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.mcool" | wc -l | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No .mcool files found in $INPUT_DIR" >&2
  exit 1
fi

echo "Submitting array: 0-$((N-1))  (samples=$N)"
sbatch --array=0-$((N-1)) "$SCRIPT"
