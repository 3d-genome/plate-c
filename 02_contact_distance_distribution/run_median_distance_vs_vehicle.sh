#!/bin/bash
#SBATCH --job-name=median
#SBATCH --partition=tttt,owners
#SBATCH --time=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

# Start from a clean module environment (critical!)
module purge
module load system
module load python/3.9.0
module load py-numpy/1.24.2_py39

cd ~/research/plateC

# -------------------------
# BASE as input argument
# Usage:
#   sbatch median.sbatch aux_data/experiment-01_plate-c_human_hek293_drug-panel-1a_24h
# Or:
#   sbatch --export=BASE=aux_data/... median.sbatch
# -------------------------
BASE="${1:-${BASE:-}}"
[[ -n "$BASE" ]] || { echo "[ERROR] BASE not provided. Usage: sbatch $0 <BASE>"; exit 1; }
[[ -d "$BASE" ]] || { echo "[ERROR] BASE directory not found: $BASE"; exit 1; }

# -------------------------
# Auto-detect rename map + vehicle file
# -------------------------
shopt -s nullglob
rename_files=( "$BASE"/rename*.txt )
veh_files=( "$BASE"/*treatment_vehicle.txt )
shopt -u nullglob

if (( ${#rename_files[@]} != 1 )); then
  echo "[ERROR] Expected exactly 1 rename*.txt in $BASE, found ${#rename_files[@]}"
  printf '  %s\n' "${rename_files[@]:-NONE}"
  exit 1
fi

if (( ${#veh_files[@]} != 1 )); then
  echo "[ERROR] Expected exactly 1 *treatment_vehicle.txt in $BASE, found ${#veh_files[@]}"
  printf '  %s\n' "${veh_files[@]:-NONE}"
  exit 1
fi

RENAME_MAP="${rename_files[0]}"
VEH="${veh_files[0]}"

echo "[INFO] BASE       = $BASE"
echo "[INFO] RENAME_MAP = $RENAME_MAP"
echo "[INFO] VEH        = $VEH"

python3 group_median_distance.py \
  --base "${BASE}" \
  --rename_map "${RENAME_MAP}" \
  --vehicle_file "${VEH}" \
  --processed_dir processed \
  --min_rep 2 \
  --out_tsv "${BASE}/median_distance_vs_vehicle_BH.tsv"

echo "[DONE] median distance finished."

