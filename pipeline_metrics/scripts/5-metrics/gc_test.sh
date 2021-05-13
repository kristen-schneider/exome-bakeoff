#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=gc_bias_stats_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/scripts/5-metrics/gc_bias_stats_driver.out
#SBATCH --error=/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/scripts/5-metrics/gc_bias_stats_driver.err

# PURPOSE: wrapping michael's script across all bam files
# michael's scirpt: /Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/scripts/5-metrics/gc_bias_stats.py

# directories
in_dir='/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/'
out_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/analyses/gcbias/59_ACMG/metric_files/'
59_ACMG_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/regions/59_ACMG/'
# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for tech_bam in `ls $in_dir`
do
    for gene_bed in `ls $59_ACMG_dir`
        do
            # bam, bed, out
            python gc_bias_stats.py $tech_bam $gene_bed $out_f
        done
done




