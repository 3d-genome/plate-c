#!/bin/bash
#
#SBATCH --job-name=map
#SBATCH --partition=tttt,owners
#SBATCH --time=0-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

fastq_r1_file_list=$1
output_parent_folder="processed"
genome_file="/home/users/tttt/data/references/hs37d5/bwa_index/genome.fa"
#genome_file="/home/users/tttt/data/references/GRCh38/bwa_index/genome.fa"
bwa_path="/home/users/tttt/tools/bwa/bwa-0.7.17"

input=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${fastq_r1_file_list})

id=${input##*/}
id=${id%%.R[12]*}
sm=${id}
echo "ID:"$id
file1=$input
file2=${input/.R1./.R2.}
echo "R1:"$file1
echo "R2:"$file2

output_folder=${output_parent_folder}/$id

mkdir ${output_folder}
ln -s $(readlink -m ${file1}) ${output_folder}/R1.fq.gz
ln -s $(readlink -m ${file2}) ${output_folder}/R2.fq.gz

${bwa_path}/bwa mem -5SP -t2 -R '@RG\tID:'"${id}"'\tPL:ILLUMINA\tSM:'"${sm}" ${genome_file} ${output_folder}/R1.fq.gz ${output_folder}/R2.fq.gz 2> ${output_folder}/aln.sam.log | gzip > ${output_folder}/aln.sam.gz


