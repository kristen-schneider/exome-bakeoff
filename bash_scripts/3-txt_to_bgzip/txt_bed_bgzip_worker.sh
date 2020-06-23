#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=txt_bed_bgzip_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/3-txt_to_bgzip/conversion_logs/txt_bed_bgzip_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/3-txt_to_bgzip/conversion_logs/txt_bed_bgzip_worker.err

# load samtools to make bgzip or tabix available on the executing node
module load samtools

echo "input file: $1"
echo "output file: $2"

bgzip -c $1 > $2
#cat $1 | awk '{OFS="\t"; print $1,$2,$2+1,$3,$4,$5,$6,$7;}' | bgzip -c > $2
