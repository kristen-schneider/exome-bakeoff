#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=cram_to_bam_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/sort_bed_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/sort_bed_driver.err

# PURPOSE: makes bam files from cram files

# directories
in_dir='/Users/krsc0813/exome-bakeoff/Analyses/quality/59_ACMG/downsampled/metric_files/'
out_dir='/Users/krsc0813/exome-bakeoff/Analyses/quality/59_ACMG/downsampled/metric_files/'

# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for f_txt in `ls $in_dir`
do
    if [[ $f_txt == *.txt ]]; then
        echo $f_txt

        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${f_txt%%.*}
        out_f=${out_dir}${sample}_sorted.txt

        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"
        else
            echo "submitting a job"
            qsub sort_bed_worker.sh $in_dir$f_txt $out_f
        fi
    fi
done
