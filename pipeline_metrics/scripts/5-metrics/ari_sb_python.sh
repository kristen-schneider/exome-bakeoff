# Ari's Attempt

#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=strandbias_python
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/armc4884/exome-bakeoff/bash_scripts/5-metrics/strandbias_python.out
#SBATCH --error=/Users/armc4884/exome-bakeoff/bash_scripts/5-metrics/strandbias_python.err
# In Quality these were put into the bash_scripts section instead of pipeline_metrics. What is the point of this?

#commandline arguments

# 1. Sample Name
sample=$1

# 2. Pileup File
pileup=$2

# 3. Regions
regions=$3

# 4. Output Directory
out_dir=$4

python strandbias.py $sample $pileup $regions $out_dir > $out_dir$sample"_strandbias.txt"

#echo $out_dir$sample"_strandbias.txt" 


