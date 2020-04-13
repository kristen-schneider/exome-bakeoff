#!/bin/bash    
#SBATCH -p long             	# Partition or queue. In this case, short!
#SBATCH --job-name=mpileup   	# Job name
#SBATCH --mail-type=NONE               	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                   	# Only use a single node
#SBATCH --ntasks=1                    	# Run on a single CPU
#SBATCH --mem=30GB                   	# Memory limit
#SBATCH --time=120:00:00               	# Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/mpileup.out 	# Standard output and error log
#SBATCH --error=/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/mpileup.err 	# %j inserts job number



module load samtools

reference='/scratch/Shares/layer/ref/hs37-1kg/human_g1k_v37.fasta'

for bam in /scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/bam_output/*.bam; do samtools mpileup -f $reference -s $bam > $bam.txt; done

