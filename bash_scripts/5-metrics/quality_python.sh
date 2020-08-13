#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=quality_python
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/quality_python.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/quality_python.err

#commandline arguments

# 1. sample name
sample=$1
# 2. pileup file
pileup=$2
# 3. regions
regions=$3
# 4. output
out_dir=$4


# make noise metric files
python quality.py $pileup $regions $out_dir > $out_dir$sample"_quality.txt"
