#!/bin/bash
#
#SBATCH --job-name=dip-c_prep
#SBATCH --partition=tttt,owners
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

dipc_path="/home/users/tttt/tools/dip-c"

folder_list=$1

folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${folder_list})

${dipc_path}/scripts/hickit_pairs_to_con.sh ${folder}/contacts_unisex.pairs.gz

