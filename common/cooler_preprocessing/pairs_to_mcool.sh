#!/usr/bin/bash
#SBATCH --job-name=pairs2mcool
#SBATCH --output=logs/pairs2mcool.%A_%a.out
#SBATCH --error=logs/pairs2mcool.%A_%a.err
#SBATCH --time=4:00:00
#SBATCH -p tttt,owners
#SBATCH --cpus-per-task=4
#SBATCH --mem=75GB

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch $0 <filename.pairs.gz>" >&2
  exit 1
fi

PAIRS_FILE="$1"

module load gcc/10.3.0
module load python/3.9.0
module load hdf5
module load zlib
module load xz
source "$HOME/tools/envs/hi_c2/bin/activate"



PROJECT="$HOME/research/plate_c/FINAL"
INPUT_DIR="${PROJECT}/data/merged/pairs_input"
OUTPUT_DIR="${PROJECT}/data/merged/input"
CHROMSIZES="${PROJECT}/genome_files/hg19.chrom.sizes"

BASE_RES=5000
RESOLUTIONS="5000,10000,25000,50000,100000,250000,500000,1000000"

mkdir -p "$OUTPUT_DIR"

input_path="${INPUT_DIR}/${PAIRS_FILE}"
sample_name="${PAIRS_FILE%.pairs.gz}"
cool_path="${OUTPUT_DIR}/${sample_name}.${BASE_RES}.cool"
cool_filtered="${OUTPUT_DIR}/${sample_name}.${BASE_RES}.filtered.cool"
mcool_path="${OUTPUT_DIR}/${sample_name}.mcool"

if [[ ! -f "$input_path" ]]; then
  echo "ERROR: File not found: $input_path" >&2
  exit 1
fi

echo "===================================="
echo "Converting: $PAIRS_FILE"
echo "Output: ${sample_name}.mcool"
echo "===================================="

# Step 1: pairs.gz -> cool at base resolution
echo "Creating base cool at ${BASE_RES}bp resolution..."
cooler cload pairs \
  -c1 2 -p1 3 -c2 4 -p2 5 \
  --zero-based \
  --assembly hg19 \
  "${CHROMSIZES}:${BASE_RES}" \
  "$input_path" \
  "$cool_path"

# Step 2: Filter to canonical autosomes only
echo "Filtering to canonical autosomes..."
python3 "${PROJECT}/scripts/merged/create_cooler.py" \
  "${cool_path}" \
  "${cool_filtered}"

# Step 3: cool -> mcool with all resolutions
echo "Zoomifying to mcool..."
cooler zoomify \
  --resolutions "$RESOLUTIONS" \
  --balance \
  --nproc "${SLURM_CPUS_PER_TASK:-4}" \
  -o "$mcool_path" \
  "${cool_filtered}"

# Clean up intermediate cool files
rm -f "$cool_path" "$cool_filtered"

echo "Successfully converted $PAIRS_FILE"
echo "Output: $mcool_path"