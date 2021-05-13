#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=noise_python
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise_python.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise_python.err

#commandline arguments

# 1. sample name
sample=$1
# 2. pileup file
pileup=$2
# 3. vcf_files
vcfs=$3
# 4. gene_bed
gene_bed=$4
# 5. output
out_dir=$5


# make noise metric files
python noise.py $sample $pileup $vcfs $gene_bed $out_dir > $out_dir$sample"_noise.txt"
