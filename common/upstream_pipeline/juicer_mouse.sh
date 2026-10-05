#!/bin/bash
#
#SBATCH --job-name=juicer_mouse
#SBATCH --partition=tttt,owners
#SBATCH --time=2-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G

juice_jar="/home/users/tttt/tools/juicer_tools/juicer_tools_1.22.01.jar"
input_file="$1"

ml system
ml biology
module load java
module load x11


java -Xmx64g -jar ${juice_jar} pre ${input_file} ${input_file/.pairs.gz/.hic} mm10
java -Xmx64g -jar "${juice_jar}" addNorm "${input_file/.pairs.gz/.hic}"
