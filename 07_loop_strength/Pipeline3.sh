#!/usr/bin/bash
#SBATCH --job-name=Pipeline3_P2LL
#SBATCH --output=logs/pipeline3_p2ll_%j.out
#SBATCH --error=logs/pipeline3_p2ll_%j.err
#SBATCH --time=04:00:00
#SBATCH -p owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB

set -euo pipefail

# 1. Load Modules
module load gcc/10.3.0
module load python/3.9.0
module load hdf5
module load zlib
module load xz

# 2. Activate Environment
source "$HOME/tools/envs/hi_c2/bin/activate"

# 3. Define Paths
PROJECT="$HOME/research/plate_c/FINAL"
INPUT_DIR="${PROJECT}/data/merged/input"              
LOOP_DIR="${PROJECT}/outputs/global_calling/Loop"     
OUTPUT_CSV="${LOOP_DIR}/Global_P2LL_Summary.csv"      
SCRIPT="${PROJECT}/scripts/merged/Calculate_P2LL.py"  

# 4. Run the P2LL Calculator
echo "Starting P2LL Quantification..."
echo "Input:  $INPUT_DIR"
echo "Loops:  $LOOP_DIR"
echo "Output: $OUTPUT_CSV"

python3 -u "$SCRIPT" \
  "$INPUT_DIR" \
  "$LOOP_DIR" \
  "hg19" \
  "$OUTPUT_CSV"

echo "Pipeline 3 Finished successfully at $(date)"