#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=strandbias_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/strandbias_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/strandbias_driver.err

# PURPOSE: Computes strandbias metric

# Commandline arguments
# 1. Pileup Files
pileups=$1

# 2. Regions
regions=$2

# 3. Output
out_dir=$3
#

echo "Starting strandbias calculation"

for pileup_file in `ls $pileups`
do
    if [[ $pileup_file == *.bed* ]] && [[ $pileup_file != *.tbi ]]; then
	echo $pileup_file

	# Prepares the name of the outpuyt file by removing the pattern '.*'
	sample=${pileup_file%%.*}
	out_file=${out_dir}${sample}
	
	# Check if the output file already exists, if so, skip the conversion
	if test -e $out_file; then
		echo "$out_file already exists, so we skip processing it"
	else
		echo "Submitting a job for" $sample
		sbatch ari_sb_worker.sh $sample $pileups$pileup_file $regions $out_dir
	fi
    #else
	#echo "not using" $pileup_file
    fi
done
