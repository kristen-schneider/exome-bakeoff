#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=cram_to_bam_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/1-cram_to_bam/conversion_logs/cram_to_bam_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/1-cram_to_bam/conversion_logs/cram_to_bam_driver.err

# PURPOSE: makes bam files from cram files

# directories
in_dir='/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/'
out_dir='/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/'

# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for f_cram in `ls $in_dir`
do
    if [[ $f_cram == *.cram ]]; then
        echo $f_cram

        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${f_cram%%.*}
        out_f=${out_dir}${sample}.bam

        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"
        else
            echo "submitting a job"
            qsub cram_to_bam_worker.sh $in_dir$f_cram $out_f
        fi
    fi
done
