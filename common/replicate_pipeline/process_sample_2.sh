#!/bin/bash
#SBATCH --job-name=analyze_cool
#SBATCH --output=logs/analyze_%A_%a.out
#SBATCH --error=logs/analyze_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH -p tttt,owners
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB

set -euo pipefail

module load gcc/10.3.0 python/3.9.0 hdf5 zlib xz
source "$HOME/tools/envs/hi_c2/bin/activate"

PROJECT="$HOME/research/plate_c/FINAL"

# Check if CSV_NAME is set
if [ -z "$CSV_NAME" ]; then
    echo "ERROR: CSV_NAME environment variable not set"
    echo "This script must be submitted via submit_pipeline.sh"
    exit 1
fi

CSV="$PROJECT/scripts/replicates/lists/$CSV_NAME"

# Read sample from CSV (skip header with +1)
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$CSV")
line=${line%$'\r'}  # Strip trailing CR from CRLF line endings
IFS=',' read -r sample_name treatment experiment panel timepoint sample_path cell_type organism genome tad_loop_ref  <<< "$line"

# Clean quotes if present
sample_name=${sample_name//\"/}
sample_path=${sample_path//\"/}
cell_type=${cell_type//\"/}
organism=${organism//\"/}
genome=${genome//\"/}
TAD_path="$PROJECT/outputs/global_calling/TAD/${tad_loop_ref//\"/}"
loop_path="$PROJECT/outputs/global_calling/Loop/${tad_loop_ref//\"/}"

echo "▶ processing $sample_name"

INDIR="$PROJECT/data/replicates/balanced/$sample_name"
INPUT_COOL="$INDIR/${sample_name}_noMT_noY_noX_5000.cool"

# Check if input exists
[[ ! -f "$INPUT_COOL" ]] && echo "ERROR: input cool file not found: $INPUT_COOL" && exit 1

# Create multi-resolution mcool
MCOOL="$INDIR/${sample_name}_noMT_noY_noX.mcool"
echo "[1] creating multi-resolution mcool"
cooler zoomify -r 10000,25000,50000,100000 "$INPUT_COOL" -o "$MCOOL"

# Balance each resolution
RESOLUTIONS=(10000 25000 50000 100000)
echo "[2] balancing resolutions"

for RES in "${RESOLUTIONS[@]}"; do
    COOL_PATH="$MCOOL::resolutions/${RES}"
    RES_KB=$((RES / 1000))
    echo "── ${RES_KB} bp ──"

    if cooler dump -H bins weight "$COOL_PATH" &>/dev/null; then
        echo "weights present – skip"
    else
        echo "balancing …"
        cooler balance "$COOL_PATH" --force
        echo "balancing done"
    fi
done

# 3) Pileup analysis
echo "[3] running pileup analysis"
cell_type_clean=$(echo "$cell_type" | tr '-' '_')
ANALYSIS_OUTDIR="$PROJECT/outputs/replicate_outputs/$cell_type_clean/$sample_name"
mkdir -p "$ANALYSIS_OUTDIR"

python3 "$PROJECT/scripts/replicates/make_pup_rep.py" \
        "$MCOOL" \
        "$cell_type" \
        "$genome" \
        "$TAD_path" \
        "$loop_path" \
        "$ANALYSIS_OUTDIR" \
        "$sample_name"

echo "✓ $sample_name done"

# Print resource usage stats
echo ""
echo "=== Resource Usage ==="
sacct -j $SLURM_JOB_ID --format=JobID,Elapsed,MaxRSS,MaxVMSize,State


