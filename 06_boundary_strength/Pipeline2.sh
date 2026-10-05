#!/usr/bin/bash
#SBATCH --job-name=Pipeline2_merged
#SBATCH --output=logs/pipeline2.%A_%a.out
#SBATCH --error=logs/pipeline2.%A_%a.err
#SBATCH --time=3:00:00
#SBATCH -p tttt,owners
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --array=0-4

set -euo pipefail

module load gcc/10.3.0
module load python/3.9.0
module load hdf5
module load zlib
module load xz
source "$HOME/tools/envs/hi_c2/bin/activate"


PROJECT="$HOME/research/plate_c/FINAL"
INPUT_DIR="${PROJECT}/data/merged/input" #I'm gonna balance them here
H5_DIR="${PROJECT}/data/merged/h5" #I'm gonna balance them here

RESOLUTIONS=(5000 10000 25000 50000 100000)
RES=${RESOLUTIONS[$SLURM_ARRAY_TASK_ID]}
RES_KB=$((RES / 1000))

TAD_BASE="${PROJECT}/outputs/global_calling/TAD"
Loop_BASE="${PROJECT}/outputs/global_calling/Loop"


# Discover samples from mcool names
mapfile -t SAMPLES < <(
  find "$INPUT_DIR" -maxdepth 1 -type f -name "*.mcool" \
    -printf "%f\n" | sed 's/\.mcool$//' | sort
)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: No .mcool files found in $INPUT_DIR" >&2
  exit 1
fi
echo "Resolution: ${RES_KB}k (task ${SLURM_ARRAY_TASK_ID})"
echo "Found ${#SAMPLES[@]} samples."


#Iterating through the samples and doing TAD calls
for sample in "${SAMPLES[@]}"; do
  MCOOL_PATH="${INPUT_DIR}/${sample}.mcool"
  H5_PATH="${H5_DIR}/${sample}_${RES_KB}k.h5"

  if [[ ! -f "$H5_PATH" ]]; then
    echo "ERROR: missing H5 for sample=$sample res=${RES_KB}k: $H5_PATH" >&2
    exit 1
  fi

  # Organize outputs per sample/resolution
  OUTDIR="${TAD_BASE}/${sample}/${RES_KB}kb"
  OUTDIR2="${Loop_BASE}/${sample}/${RES_KB}kb"

  mkdir -p "$OUTDIR"
  mkdir -p "$OUTDIR2"

  TAG="${sample}_${RES_KB}kb_min$((RES*3/1000))kb_max$((RES*5/1000))kb_step$((RES/1000))kb_thr001_fdr001"
  OUTPREFIX="${OUTDIR}/TAD_${TAG}"
  OUTPREFIX2="${OUTDIR2}/Loop_${TAG}"

  echo "-------------------------------------"
  echo "Sample: $sample | RES: ${RES_KB}k"
  echo "H5: $H5_PATH"
  echo "OUTPREFIX: $OUTPREFIX"
  echo "-------------------------------------"

  echo "[TAD] hicFindTADs"
  hicFindTADs \
    -m "$H5_PATH" \
    --outPrefix "$OUTPREFIX" \
    --thresholdComparisons 0.01 \
    --delta 0.01 \
    --minDepth $((RES * 3)) \
    --maxDepth $((RES * 5)) \
    --step "$RES" \
    --correctForMultipleTesting fdr

  echo "[Loop] Loop_finder.py"
  python3 -u "$PROJECT/scripts/merged/Loop_finder.py" \
    "$MCOOL_PATH" "$OUTDIR" "$OUTPREFIX" "$RES" "$sample"

  echo "Done sample=$sample res=${RES_KB}k at $(date)"
done

echo "Finished array task ${SLURM_ARRAY_TASK_ID} (RES=${RES_KB}k) at $(date)"
