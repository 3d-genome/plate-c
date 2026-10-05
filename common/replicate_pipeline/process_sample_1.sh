#!/bin/bash
#SBATCH --job-name=make_cool
#SBATCH --output=logs/cool_%A_%a.out
#SBATCH --error=logs/cool_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH -p tttt,owners
#SBATCH --cpus-per-task=4
#SBATCH --mem=75G

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

# Read sample from CSV
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$CSV")
line=${line%$'\r'}  # Strip trailing CR from CRLF line endings
IFS=',' read -r sample_name treatment experiment panel timepoint sample_path cell_type organism genome tad_loop_ref <<< "$line"

# Clean quotes if present
sample_name=${sample_name//\"/}
sample_path=${sample_path//\"/}
cell_type=${cell_type//\"/}
organism=${organism//\"/}
genome=${genome//\"/}

echo "▶ processing $sample_name"

# Genome setup
if [[ "$organism" == "human" ]]; then
    GENOME_FILE="$PROJECT/genome_files/hg19.chrom.sizes"
    GENOME_BUILD="hg19"
elif [[ "$organism" == "mouse" ]]; then
    GENOME_FILE="$PROJECT/genome_files/mm10.chrom.sizes"
    GENOME_BUILD="mm10"
else
    echo "ERROR: Unknown organism: $organism"; exit 1
fi

# Paths
PAIRS="${sample_path}/contacts_unisex.pairs.gz"
[[ ! -f "$PAIRS" ]] && echo "ERROR: pairs not found: $PAIRS" && exit 1

# Extract sample folder name from path (everything after 'processed/')
OUTDIR="$PROJECT/data/replicates/balanced/$sample_name"
mkdir -p "$OUTDIR"

TEMP_COOL="$OUTDIR/${sample_name}_temp.cool"
FINAL_COOL="$OUTDIR/${sample_name}_noMT_noY_noX_5000.cool"

# Skip if done
[[ -f "$FINAL_COOL" ]] && echo "✓ already exists, skipping" && exit 0

# ==== NEW PART: pairs.gz → cool ====
# [1] Step removed (no need to make manual bins)
echo "[2] pairs → cool"
# Use GENOME_FILE directly with the :5000 suffix
cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    --zero-based \
    --assembly "$GENOME_BUILD" \
    "$GENOME_FILE:5000" \
    "$PAIRS" \
    "$TEMP_COOL"

# echo "[1] creating bins (5kb)"
# BINS="$OUTDIR/${sample_name}_bins.bed"
# cooler makebins "$GENOME_FILE" 5000 > "$BINS"
# echo "[2] pairs → cool"
# cooler cload pairs \
#     -c1 2 -p1 3 -c2 4 -p2 5 \
#     --assembly "$GENOME_BUILD" \
#     "$BINS:5000" \
#     "$PAIRS" \
#     "$TEMP_COOL"

# ==== YOUR EXISTING CODE: filter chromosomes ====
echo "[3] filtering chromosomes"
python3 "$PROJECT/scripts/replicates/create_cooler.py" "$TEMP_COOL" "$FINAL_COOL"

# Cleanup
rm -f "$TEMP_COOL" #"$BINS"

echo "✓ done: $FINAL_COOL"
ls -lh "$FINAL_COOL"

# Print resource usage stats
echo ""
echo "=== Resource Usage ==="
sacct -j $SLURM_JOB_ID --format=JobID,Elapsed,MaxRSS,MaxVMSize,State


