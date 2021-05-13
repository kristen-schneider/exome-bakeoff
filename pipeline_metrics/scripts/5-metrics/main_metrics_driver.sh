#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=main_metric_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/main_metrics_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/main_metrics_driver.err


# commandline arguments
# 1. metric
metric=$1

# quality
if [[ $metric == 'quality' ]]; then
    pileups=$2
    region=$3
    out_dir=$4
    echo "computing" $metric
    sbatch quality_driver.sh $pileups $region $out_dir
fi

# strand bias
if [[ $metric == 'sb' ]]; then
    echo $metric
fi

# noise
if [[ $metric == 'noise' ]]; then
    pileups=$2
    vcfs=$3
    region=$4
    out_dir=$5
    echo "computing" $metric
    sbatch noise_driver.sh $pileups $vcfs $region $out_dir
fi
