#!/bin/bash
#SBATCH --job-name=dipc_color2
#SBATCH --partition=tttt,owners
#SBATCH --time=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

dipc_path="/home/users/tttt/tools/dip-c"
infile="$1"
outfile="${infile%.con.gz}.cpg_b1m.color2"

ml python/2.7.13
ml py-numpy/1.14.3_py27
ml py-scipy/1.1.0_py27

echo "[$(date)] color2: $infile -> $outfile"
${dipc_path}/dip-c color2 -b1000000 -H -c ${dipc_path}/color/hg19.cpg.1m.txt -s "$infile" > "$outfile"

