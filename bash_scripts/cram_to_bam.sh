#!/usr/bin/env bash
#
#SBATCH -p long
#SBATCH --job-name=cram_to_bam
#SBATCH --ntasks=1
#SBATCH --time=40:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=cram_to_bam.out
#SBATCH --error=cram_to_bam.err

module load samtools

ref=/Shares/layer_shared/ref/hg37/human_g1k_v37.fasta.gz
in_cram_dir=/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/
#in_cram=/Users/krsc0813/exome-bakeoff/AgilentQXT-IDT-0720ME25_S12_L001.cram
out_bam_dir=/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/
#out_bam=/Users/krsc0813/exome-bakeoff/AgilentQXT-IDT-0720ME25_S12_L001.bam

for cram in `ls ${in_cram_dir}*.cram`
do

    base_name=`basename $cram`
    sample=${base_name%%.cram}
    #echo $out_bam_dir${sample}.cram
    samtools view -T $ref -b -o $out_bam_dir${sample}.bam $in_cram_dir${sample}.cram \
        && samtools index $out_bam_dir${sample}.bam
done
