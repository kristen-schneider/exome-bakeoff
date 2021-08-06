#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=strandbias_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/armc4884/exome-bakeoff/bash_scripts/5-metrics/strandbias_worker.out 
#SBATCH --error=/Users/armc4884/exome-bakeoff/bash_scripts/5-metrics/strandbias_worker.err

#module load samtools
#module load pysam

# Commandline arguemnts
# 1. Sample Name
sample=$1

# 2. Pileup File
pileup=$2

# 3. Regions
regions=$3

# 4. Output
out_dir=$4

sbatch strandbias_python.sh $sample $pileup $regions $out_dir
