#!/bin/bash
#
#SBATCH --job-name=juicer_human
#SBATCH --partition=tttt,owners
#SBATCH --time=2-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=128G

juice_jar="/home/users/tttt/tools/juicer_tools/juicer_tools_1.22.01.jar"
input_file="$1"
ml system
ml biology
module load java
module load x11


java -Xmx16g -jar ${juice_jar} pre ${input_file} ${input_file/.pairs.gz/.hic} hg19
java -Xmx16g -jar "${juice_jar}" addNorm "${input_file/.pairs.gz/.hic}"

