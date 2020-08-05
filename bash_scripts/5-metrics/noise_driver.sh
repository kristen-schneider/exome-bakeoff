#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=noise_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/noise_driver.err

# PURPOSE: computes noise metric

# commandline arguments
# 1. pileup_files
pileups=$1
# 2. vcf_files
vcfs=$2
# 3. regions
regions=$3
# 4. output
out_dir=$4
#in_dir='/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/'
#out_dir='/Users/krsc0813/exome-bakeoff/Analyses/noise/59_ACMG'
#
echo "starting noise calculation"

for pileup_file in `ls $pileups`
do
    if [[ $pileup_file == *.bed* ]] && [[ $pileup_file != *.tbi ]]; then
        echo $pileup_file
        
        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${pileup_file%%.*}
        out_file=${out_dir}${sample}
       
        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_file; then
            echo "$out_file already exists, so we skip processing it"
        else
            echo "submitting a job for" $sample
            qsub noise_worker.sh $sample $vcfs $regions $out_dir
        fi 
    #else
        #echo "not using" $pileup_file
    fi
done
