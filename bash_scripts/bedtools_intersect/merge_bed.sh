#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=merge_bed
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/merge.out
#SBATCH --error=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/merge.err

# PURPOSE: merge bed file 

module load bedtools

output_intersection_bed=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/output_intersect_sorted.bed
#output_intersection_bed=/Shares/layer_shared/projects/chco/exome_bakeoff/regions/intersect_bed_vcf/test.bed

bedtools merge -i $output_intersection_bed
