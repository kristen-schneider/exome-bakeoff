#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=bam_to_mpileup_driver
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=bam_to_mpileup_driver.out
#SBATCH --error=bam_to_mpileup_driver.err

# Iterate through the dir
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
            echo "$out_f already exists, so skipping processin it"
        else
            echo "submitting a job"
            qsub bam_to_mpileup_worker.sh $in_dir$f_bam $out_f
        fi
    fi
done


