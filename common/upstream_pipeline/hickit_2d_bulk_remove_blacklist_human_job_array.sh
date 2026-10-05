#!/bin/bash
#
#SBATCH --job-name=hickit_2d
#SBATCH --partition=tttt,owners
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

PATH="/home/users/tttt/bin:$PATH"
hickit_path="/home/users/tttt/tools/hickit/"
blacklist_file="hg19-blacklist.v2-renamed.bed"

folder_list=$1

folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${folder_list})


## create segment files
# unisex
${hickit_path}/hickit.js sam2seg ${folder}/aln.sam.gz 2> ${folder}/contacts_unisex.seg.log | ${hickit_path}/hickit.js chronly - | ${hickit_path}/hickit.js bedflt ${blacklist_file} - | gzip > ${folder}/contacts_unisex.seg.gz

## creat contact files
# unisex
${hickit_path}/hickit --dup-dist=1 -i ${folder}/contacts_unisex.seg.gz -o - 2> ${folder}/contacts_unisex.pairs.log | bgzip > ${folder}/contacts_unisex.pairs.gz

## summarize information
sample=${folder##*/} # sample name
dup_line=$(grep "duplicate" ${folder}/contacts_unisex.pairs.log)
dup=${dup_line%%\%*};dup=${dup##* } # dup rate
dup_num=${dup_line%% /*};dup_num=${dup_num##* }
raw=${dup_line##* } # raw contacts
con=$((raw-dup_num)) # contacts
intra=$(zcat ${folder}/contacts_unisex.pairs.gz | grep -v "^#" | awk '{sum++;if($2==$4){intra++}}END{print intra*100/sum}') # percent intra
echo ${sample} ${raw} ${dup} ${con} ${intra} > ${folder}/contacts_unisex.info


