
#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=strandbias_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/strandbias_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/strandbias_driver.err

# PURPOSE: driver script to iterate over directories of mpileup/regions files

# commandline arguments
# 1. pileup_files
pileups=$1
# 2. out_dir
out_dir=$2

echo "starting strandbias calculations..."

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
            sbatch quality_worker.sh $sample $pileups$pileup_file $regions $out_dir
        




