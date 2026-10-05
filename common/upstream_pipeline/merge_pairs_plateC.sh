#!/bin/bash
#
#SBATCH --job-name=merge_pairs
#SBATCH --partition=tttt,owners
#SBATCH --time=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

header_file="$1"
cellname_file="$2"
out_file="$3"

PATH="/home/users/tttt/bin:$PATH"

echo "generating: ${out_file}.gz"
cat ${header_file} > ${out_file}
for f in `cat ${cellname_file} | awk '{print $0"/contacts_unisex.pairs.gz"}'`; do >&2 echo "  appending: $f"; gunzip -c $f | grep -v "^#" | cut -f1-7; done >> ${out_file}
echo "  merging and compressing"
hickit --dup-dist=0 -i ${out_file} -o - | bgzip > ${out_file}.gz
rm ${out_file}

