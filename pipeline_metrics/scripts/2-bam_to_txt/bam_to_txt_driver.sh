#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=bam_to_txt_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/2-bam_to_txt/conversion_logs/bam_to_txt_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/2-bam_to_txt/conversion_logs/bam_to_txt_driver.err

# PURPOSE: makes mpileup files from bam files

# directories
in_dir='/Shares/layer_shared/projects/sequence_analysis/kristen_downsample/align/'
out_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_downsample/'

# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for f_bam in `ls $in_dir` 
do
    if [[ $f_bam == *.bam ]]; then
        echo $f_bam

        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${f_bam%%.*}
        out_f=${out_dir}${sample}.txt
        
        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"
        else
            echo "submitting a job"
            qsub bam_to_txt_worker.sh $in_dir$f_bam $out_f
        fi
    fi
done


