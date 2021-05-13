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

set -e

# directories
in_dir=/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/
out_dir=/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/analyses/gcbias/oncogenes/metric_files/
region_file=/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/regions/all_oncogenes.bed
ref=/Shares/layer_shared/ref/hs37-1kg/human_g1k_v37.fasta
# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for tech_bam in `ls $in_dir`
do
    if [[ $tech_bam == *.bam ]]; then
        echo $tech_bam 
        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${tech_bam%%.*}
        out_f=${out_dir}${sample}_gcbias.txt

        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"
            # bam, bed, out
        else
            python gc_bias_stats.py $in_dir$tech_bam $region_file $out_f $ref
            #echo "sbatch gc_bias_stats_worker.sh $in_dir$tech_bam $region_dir$gene_bed $out_f"
            #python gc_bias_stats.py $in_dir$tech_bam $region_dir$gene_bed $out_f
    
        fi
    fi
done




