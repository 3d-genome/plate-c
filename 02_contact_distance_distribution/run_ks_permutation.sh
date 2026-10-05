#!/bin/bash
#SBATCH --job-name=groupKS
#SBATCH --partition=tttt,owners
#SBATCH --time=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

module purge
module load system
module load python/3.9.0
module load py-numpy/1.24.2_py39
# module load py-scipy/<..._py39>

cd ~/research/plateC

# -------------------------
# BASE as input argument
# Usage:
#   sbatch groupKS.sbatch aux_data/experiment-01_plate-c_human_hek293_drug-panel-1a_24h
# Or:
#   sbatch --export=BASE=aux_data/... groupKS.sbatch
# -------------------------
BASE="${1:-${BASE:-}}"
[[ -n "$BASE" ]] || { echo "[ERROR] BASE not provided. Usage: sbatch $0 <BASE>"; exit 1; }
[[ -d "$BASE" ]] || { echo "[ERROR] BASE directory not found: $BASE"; exit 1; }

# -------------------------
# Auto-detect RENAME_MAP and VEH
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

OUTDIR="${BASE}/groupKS"
mkdir -p "${OUTDIR}"

echo "[INFO] BASE       = $BASE"
echo "[INFO] RENAME_MAP = $RENAME_MAP"
echo "[INFO] VEH        = $VEH"
echo "[INFO] OUTDIR     = $OUTDIR"

K_PER_REP=200000
GRID_N=2000
N_PERM=20000
SEED=1
MIN_REP=2
MAX_FULL_PERM=200000

# -------------------------
# 1) group-level KS for each treatment (skip vehicle)
# -------------------------
shopt -s nullglob
#treat_files=( "$BASE"/*_treatment_*.txt )
treat_files=( "$BASE"/experiment*_cluster_*.txt )
shopt -u nullglob

if (( ${#treat_files[@]} == 0 )); then
  echo "[ERROR] No treatment files found in $BASE matching *_treatment_*.txt"
  exit 1
fi

for f in "${treat_files[@]}"; do
  [[ "$f" == "$VEH" ]] && continue

  bn="$(basename "$f")"
  tag="${bn%.txt}"
  out="${OUTDIR}/${tag}_vs_vehicle.tsv"

  echo "[INFO] ${bn} vs vehicle -> ${out}"

  python3 group_ks_ecdf_perm_with_sign.py \
    --groupA_file "$f" \
    --groupB_file "$VEH" \
    --rename_map "$RENAME_MAP" \
    --processed_dir processed \
    --out_tsv "$out" \
    --k_per_rep "${K_PER_REP}" \
    --grid_n "${GRID_N}" \
    --n_perm "${N_PERM}" \
    --seed "${SEED}" \
    --min_rep "${MIN_REP}" \
    --max_full_perm "${MAX_FULL_PERM}"
done

# -------------------------
# 2) summarize results
# -------------------------
python3 summarize_groupKS.py \
  --groupks_dir "${OUTDIR}" \
  --out_tsv "${BASE}/groupKS_summary_vs_vehicle_with_sign.tsv"


# -------------------------
# 2) calculate BH-FDR
# -------------------------
python3 bh_fdr_groupKS.py \
  --groupks_dir "${BASE}/groupKS" \
  --out_tsv "${BASE}/groupKS_summary_vs_vehicle_BH_with_sign.tsv"

echo "[DONE] All group KS + summary finished."

