#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=intersect_regions_ss
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/splice_sites
#SBATCH --error=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/splice_sites

# PURPOSE: intersect all splice sites with all bed files from a given region

module load bedtools


#echo $1
#echo $2
#echo $3
bedtools intersect  -a $1 -b $2 > $3

#SCAP_vcf=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/scap_COMBINED_v1.0.vcf

#exons_bed=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/test.bed
#SCAP_vcf=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/test.vcf

#bedtools intersect  -a $exons_bed -b $SCAP_vcf > output_intersect.bed
