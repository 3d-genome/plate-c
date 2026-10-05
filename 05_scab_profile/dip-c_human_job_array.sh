#!/bin/bash
#
#SBATCH --job-name=dip-c
#SBATCH --partition=tttt,owners
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

dipc_path="/home/users/tttt/tools/dip-c"

folder_list=$1

folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${folder_list})

# load modules
ml python/2.7.13
ml py-numpy/1.14.3_py27
ml py-scipy/1.1.0_py27

# compartment for PCA
${dipc_path}/dip-c color2 -b1000000 -H -c ${dipc_path}/color/hg19.cpg.1m.txt -s ${folder}/contacts_unisex.con.gz > ${folder}/cpg_b1m.color2

