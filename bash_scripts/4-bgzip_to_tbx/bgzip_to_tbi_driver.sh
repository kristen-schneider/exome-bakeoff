#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=bgzip_to_tbi_driver
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/4-bgzip_to_tbx/conversion_logs/bgzip_to_tbi_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/4-bgzip_to_tbx/conversion_logs/bgzip_to_tbi_driver.err

#pip install pysam

in_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_full_tbi/'
out_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_full_tbi/'

# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for f_gz in `ls $in_dir`
do
    if [[ $f_gz == *.bed.gz ]]; then
        echo $f_gz

        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${f_gz%%.*}
        out_f=${out_dir}${sample}.tbi

        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"
        else
            echo "submitting a job"
            sbatch bgzip_to_tbi_worker.sh $in_dir$f_gz
        fi
    fi
done
