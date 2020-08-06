#!/usr/bin/env bash
#
#SBATCH -p short
#SBATCH --job-name=quality_worker
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/quality_worker.out
#SBATCH --error=/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/quality_worker.err


#module load samtools
#module load pysam

# commandline arguments
# 1. sample name
sample=$1
# 2. pileup file
pileup=$2
# 3. regions
regions=$3
# 4. output
out_dir=$4

sbatch quality_python.sh $sample $pileup $regions $out_dir

#for vcf_file in `ls $vcfs`
#do
    #echo $pileup
    #echo $vcfs
    #echo $regions
    #echo $out_dir        

#    if [[ $vcf_file == *"$sample"* ]]; then 
#        echo $sample
#        echo $pileup
#        echo $vcf_file
#        
#        sbatch noise_python.sh $sample $pileup $vcfs $regions $out_dir
#        #python noise.py $pileup $vcfs $regions $out_dir > $sample.txt 
#    fi
#done
#ref=/Shares/layer_shared/ref/hg37/human_g1k_v37.fasta.gz

#samtools view -T $ref -b -o $2 $1 \
    #&& samtools index $2
