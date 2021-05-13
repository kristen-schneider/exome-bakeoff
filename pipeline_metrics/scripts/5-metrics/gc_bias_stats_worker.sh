#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=gc_bias_stats_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/scripts/5-metrics/gc_bias_stats_worker.out
#SBATCH --error=/Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/scripts/5-metrics/gc_bias_stats_worker.err

# PURPOSE: wrapping michael's script across all bam files
# michael's scirpt: /Shares/layer_shared/projects/chco/exome_bakeoff/pipeline/scripts/5-metrics/gc_bias_stats.py

set -e
bam=$1
bed=$2
out=$3
ref=$4

# bam, bed, out
python gc_bias_stats.py $bam $bed $out $ref


