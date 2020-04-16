#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=txt_bed_bgzip_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/3-txt_to_bgzip/conversion_logs/txt_bed_bgzip_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/3-txt_to_bgzip/conversion_logs/txt_bed_bgzip_driver.err

in_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_downsample/'
out_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_downsample/'

# create the out dir if not exists
test ! -d $out_dir && mkdir -p $out_dir

for f_txt in `ls $in_dir`
do
    if [[ $f_txt == *.txt ]]; then
        echo $f_txt

        # prepare the name of the output file by removing the pattern '.*' greedily
        sample=${f_txt%%.*}
        out_f=${out_dir}${sample}.bed.gz
    
        # Check if the output file already exists, if so, skip the conversion
        if test -e $out_f; then
            echo "$out_f already exists, so skipping processing it"
        else
            echo "submitting a job"
            qsub txt_bed_bgzip_worker.sh $in_dir$f_txt $out_f
        fi
    fi
done
