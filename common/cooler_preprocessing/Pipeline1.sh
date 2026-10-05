#!/usr/bin/bash
#SBATCH --job-name=Pipeline1_merged
#SBATCH --output=logs/pipeline1.%A_%a.out
#SBATCH --error=logs/pipeline1.%A_%a.err
#SBATCH --time=4:00:00
#SBATCH -p tttt,owners
#SBATCH --cpus-per-task=4
#SBATCH --mem=75GB

set -euo pipefail

module load gcc/10.3.0
module load python/3.9.0
module load hdf5
module load zlib
module load xz
source "$HOME/tools/envs/hi_c2/bin/activate"

PROJECT="$HOME/research/plate_c/FINAL"
INPUT_DIR="${PROJECT}/data/merged/input/" #I'm gonna balance them here
OUTPUT_DIR="${PROJECT}/data/merged/h5/" #I'm gonna balance them here

# Drop Y and MT for the 5kb resolution
# echo "[0] Dropping MT and Y and X"
# python3 create_cooler.py "${PROJECT}/processed/withMT/merged_primary_HEK.mcool::resolutions/5000" "${PROJECT}/processed/noMT_noY_noX/merged_primary_HEK_noMT_noY_noX_5000.cool"
# INPUT_COOL="${PROJECT}/processed/noMT_noY_noX/merged_primary_HEK_noMT_noY_noX_5000.cool"

# --- Auto-discover samples and select one based on array task ID ---
mapfile -t SAMPLES < <(
  find "$INPUT_DIR" -maxdepth 1 -type f -name "*.mcool" \
  -printf "%f\n" | sed 's/\.mcool$//' | sort
)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: No .mcool files found in $INPUT_DIR" >&2
  exit 1
fi

echo "Found ${#SAMPLES[@]} samples total"

# Select sample based on SLURM_ARRAY_TASK_ID
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID not set. This script must be run as an array job." >&2
  exit 1
fi

sample="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
OUTPUT_MCOOL="${INPUT_DIR}/${sample}.mcool"

echo "===================================="
echo "Processing sample: $sample (array task $SLURM_ARRAY_TASK_ID)"
echo "===================================="

#Running the balancing script on these resolutions for this sample
RESOLUTIONS=(5000 10000 25000 50000 100000)

# echo "Balancing..."
# for r in "${RESOLUTIONS[@]}"; do
#   echo "  [balance] ${r}"
#   cooler balance "${OUTPUT_MCOOL}::resolutions/${r}"
# done

echo "Filtering and converting to .h5..."
for r in "${RESOLUTIONS[@]}"; do
  filtered_cool="${OUTPUT_DIR}/${sample}_$((r/1000))k_filtered.cool"
  out_h5="${OUTPUT_DIR}/${sample}_$((r/1000))k.h5"

  echo "  [filter] ${r} → canonical autosomes only"
  python3 "${PROJECT}/scripts/merged/create_cooler.py" \
    "${OUTPUT_MCOOL}::resolutions/${r}" \
    "${filtered_cool}"

  echo "  [h5] ${r} → $(basename "$out_h5")"
  hicConvertFormat \
    --matrices "${filtered_cool}" \
    --inputFormat cool \
    --outputFormat h5 \
    --outFileName "$out_h5"

  rm -f "${filtered_cool}"
done

echo "Sample $sample processed successfully."
