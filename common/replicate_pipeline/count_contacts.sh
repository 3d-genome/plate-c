#!/bin/bash
#SBATCH --job-name=count_contacts
#SBATCH --output=logs/count_%j.out
#SBATCH --error=logs/count_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH -p tttt,owners

module load python/3.9.0

python3 count_contacts.py lists/qc_passed_samples.csv