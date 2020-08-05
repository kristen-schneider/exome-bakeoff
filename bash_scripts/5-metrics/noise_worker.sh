#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=noise_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise/noise_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise/noise_worker.err


module load samtools
module load pysam

# commandline arguments
# 1. sample name
sample=$1
# 2. vcf_files
vcfs=$2
# 3. regions
regions=$3
# 4. output
out_dir=$4

for vcf_file in `ls $vcfs`
do
    if [[ vcf_file == *$sample* ]]; then 
        echo $vcf_file
    fi
done
#ref=/Shares/layer_shared/ref/hg37/human_g1k_v37.fasta.gz

#samtools view -T $ref -b -o $2 $1 \
    #&& samtools index $2
