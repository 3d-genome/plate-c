#!/bin/bash
set -euo pipefail

BASE=$1
ml python/3.9.0

python3 bh_fdr_groupKS.py \
  --groupks_dir "${BASE}/groupKS" \
  --out_tsv "${BASE}/groupKS_summary_vs_vehicle_BH.tsv"

