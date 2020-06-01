#!/bin/bash

seqtk="/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/seqtk/seqtk"
fifty_million=50000000
sub_dir="/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/subsample/"

for sample in /Shares/layer_shared/projects/chco/exome_bakeoff/fastq_files/*
do
    if [ -d "$sample" ]
    then
        dir="/*"
        sample_dir=$sample$dir
        
        for read_num in $sample_dir
        do
            read_name=$(basename -s .fastq.gz $read_num)
            extension=".fq"
            read_name_fq=$sub_dir$read_name$extension
            $seqtk sample -s100 $read_num $fifty_million > $read_name_fq
        done
    fi
done

