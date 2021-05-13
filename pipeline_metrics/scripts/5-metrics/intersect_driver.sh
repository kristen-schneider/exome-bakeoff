#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=intersect_driver
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/intersect_driver.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/intersect_driver.err

# PURPOSE: computes quality metric

# commandline arguments
# 1. tech list
pairs_dir=$1
# 2. vcf files
vcf_dir=$2
# 3. out directory
out_dir=$3

echo "starting  calculation"

for tech_pairing in `ls $pairs_dir`
do
    sbatch intersect_worker.sh $pairs_dir $tech_pairing $vcf_dir $out_dir
done


#while IFS= read -r line; do
#    echo $line
#    tech=$line
#    
#    for vcf_file in `ls $vcf_dir`
#    do
#        #if [[ $vcf_file == *.gz* ]] && [[ $vcf_file != *.tbi ]] && [[ vcf_file == $tech* ]]; then
#        if [[ $vcf_file == *$tech* ]] && [[ $vcf_file == *.gz ]]; then
#            out_file=$out_dir$tech.txt
#            
#            # Check if the output file already exists, if so, skip the conversion
#            if test -e $out_file; then
#                echo "$out_file already exists, so we skip processing it"
#            else
#                echo "submitting a job for" $tech
#                sbatch intersect_worker.sh $tech $vcf_dir $out_file
#            fi 
#        fi
#    done
#done < $tech_list

        
        # prepare the name of the output file by removing the pattern '.*' greedily
        #sample_technology=${vcf_file%%-*}
        #out_file=${out_dir}${sample}
#        echo $sample
#        #sbatch quality_worker.sh $sample $pileups$pileup_file $regions $out_dir
