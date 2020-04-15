#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=cram_to_bam_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/1-cram_to_bam/conversion_logs/cram_to_bam_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/1-cram_to_bam/conversion_logs/cram_to_bam_driver.err


module load samtools

ref=/Shares/layer_shared/ref/hg37/human_g1k_v37.fasta.gz

samtools view -T $ref -b -o $2 $1 \
    && samtools index $2
