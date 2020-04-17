#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=bgzip_to_tbi_worker
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/4-bgzip_to_tbx/conversion_logs/bgzip_to_tbi_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/4-bgzip_to_tbx/conversion_logs/bgzip_to_tbi_worker.err

echo "input file: $1"
echo "output file: $2"

tabix $1

