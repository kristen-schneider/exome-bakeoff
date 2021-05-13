#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=noise_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise_worker.err


module load samtools
#module load pysam

# commandline arguments
# 1. sample name
sample=$1
# 2. pileup file
pileup=$2
# 3. vcf_files
vcfs=$3
# 4. regions
regions=$4
# 5. output
out_dir=$5

#for gene_bed in `ls $regions`
#do
sbatch noise_python.sh $sample $pileup $vcfs $regions $out_dir
#done
