#!/bin/bash
#
#SBATCH --job-name=hist
#SBATCH --partition=tttt,owners
#SBATCH --time=2-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

# -------------------------
# Args / params
# -------------------------
folder_list=$1

min_distance=1000
binw=0.05
ld_min=2.95
ld_max=8.35

mkdir -p logs

if [[ -z "${folder_list:-}" ]]; then
  echo "ERROR: missing folder list argument"
  echo "Usage: sbatch ... $0 folder_list_experiments.txt"
  exit 1
fi

input=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${folder_list}")

if [[ -z "${input}" ]]; then
  echo "ERROR: Empty line for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

pairs_gz="${input}/contacts_unisex.pairs.gz"
out_file="${input}/contacts_unisex.distance_log10_cis_histogram.txt"

if [[ ! -f "${pairs_gz}" ]]; then
  echo "ERROR: Missing input file: ${pairs_gz}"
  exit 1
fi

echo "Job ${SLURM_JOB_ID:-NA}, task ${SLURM_ARRAY_TASK_ID}"
echo "Folder: ${input}"
echo "Input : ${pairs_gz}"
echo "Output: ${out_file}"
echo "NOTE: cis-only histogram; trans excluded from normalization"

# -------------------------
# Build normalized cis histogram
# -------------------------
gunzip -c "$pairs_gz" \
  | grep -v "^#" \
  | awk -F $'\t' \
      -v min_distance="${min_distance}" \
      -v binw="${binw}" \
      -v ld_min="${ld_min}" \
      -v ld_max="${ld_max}" '
      BEGIN { total = 0 }

      # autosomes only
      ($2!="X"&&$2!="Y"&&$2!="chrX"&&$2!="chrY"&&
       $4!="X"&&$4!="Y"&&$4!="chrX"&&$4!="chrY") {

        # cis only
        if ($2 != $4) next

        d = $5 - $3
        if (d >= min_distance) {
          ld = log(d)/log(10)
          if (ld < ld_min || ld > ld_max) next

          b = binw * int(ld / binw)
          if (b == 0) b = 0

          hist[b]++
          total++
        }
      }

      END {
        if (total == 0) exit
        for (b in hist) {
          printf "%g %.8g\n", b, hist[b]/total
        }
      }
    ' \
  | sort -n > "$out_file"

echo "Done."

