#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=bam_to_mpileup_worker
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/2-bam_to_txt/conversion_logs/bam_to_mpileup_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/2-bam_to_txt/conversion_logs/bam_to_mpileup_worker.err


# load samtools for mpileup
module load samtools

# reference 
ref='/Shares/layer_shared/ref/hg37/human_g1k_v37.fasta.gz'

# $0 is technically the zeroth arguemnt and is set to the name of the script, 'worker.sh' in this case
# $1 will be the first argument, the input txt file in this case
# $2 will be the second argument, the to be generated output file
echo "input file: $1"
echo "output file: $2"

samtools mpileup -f $ref -s -d 8000 $1 > $2
