#!/usr/bin/env bash
#
#SBATCH -p long
#SBATCH --job-name=bam_to_mpileup_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=bam_to_mpileup_worker.out
#SBATCH --error=bam_to_mpileup_worker.err


# load samtools for mpileup
module load samtools

# reference 
ref='/Shares/layer_shared/ref/hg37/human_g1k_v37.fasta.gz'

# $0 is technically the zeroth arguemnt and is set to the name of the script, 'worker.sh' in this case
# $1 will be the first argument, the input txt file in this case
# $2 will be the second argument, the to be generated output file
echo "input file: $1"
echo "output file: $2"

samtools mpileup -f $ref -s $1 > $2
