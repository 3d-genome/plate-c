#!/usr/bin/bash
set -euo pipefail

PROJECT="$HOME/research/plate_c/FINAL"
INPUT_DIR="${PROJECT}/data/merged/input"

NSAMPLES=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.mcool" | wc -l | tr -d ' ')
if [[ "$NSAMPLES" -eq 0 ]]; then
  echo "ERROR: No .mcool files found in $INPUT_DIR" >&2
  exit 1
fi

MAX=$((NSAMPLES - 1))
echo "Submitting array 0-${MAX} (${NSAMPLES} samples total)"
sbatch --array=0-"$MAX" Pipeline1.sh
