#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=intersect_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/intersect_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/intersect_worker.err


module load bedtools

# commandline arguments
# 1. pairs dir
pairs_dir=$1
# 2. tech_sample
tech_sample=$2
# 3. vcf dir
vcf_dir=$3
# 4. output file
out_dir=$4
# 4. output
#out_dir=$4

while IFS= read -r line; do
    tech_sample_name=${tech_sample%%.*}
    sample_1="$(cut -d',' -f1 <<<$line)"
    sample_2="$(cut -d',' -f2 <<<$line)"
    out_file=$out_dir$tech_sample_name.txt
    #echo "$sample_1" >> $out_dir$tech_sample_name.txt
    echo $sample_1 $sample_2 >> $out_file
    bedtools intersect -u -a $vcf_dir$sample_1 -b $vcf_dir$sample_2 | wc -l >> $out_file
done < $pairs_dir$tech_sample



#for vcf_file_1 in `ls $vcf_dir`
#do
#    if [[ $vcf_file_1 == *$tech_name* ]] && [[ $vcf_file_1 == *.gz ]]; then
#        echo "vcf1: " $vcf_file_1 >> $out_file
#        for vcf_file_2 in `ls $vcf_dir`
#        do
#            if [[ $vcf_file_2 == *$tech_name* ]] && [[ $vcf_file_2 == *.gz ]]; then
#                echo "vcf2: " $vcf_file_2 >> $out_file
#                echo 'stop' >> $out_file
#                #bedtools intersect -u -a $vcf_dir$vcf_file_1 -b $vcf_dir$vcf_file_2 | wc -l >> $out_file
#            fi
#        done
#    fi
#done
