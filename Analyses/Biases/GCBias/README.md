# acmg_all_59.bed  
This bed file is a combination of all bed files for the acmg 59 it is created by running make_acmg_all_59_bed.sh

# filter_bams_by_acmg_59.sh
This script loops over all the bam files passing them into filter_single_bam.sh one at a time to be filtered to only contain the information inside the acmg 59 genomic ranges

# filtering_job.sh  
This is a slurm script that starts filter_bams_by_acmg_59.sh

# filter_single_bam.sh
This script takes a single bam file and filters it using samtools and acmg_all_59.bed

# make_acmg_all_59_bed.sh
This file concatonates all of the acmg 59 bed files into one. Produces acmg_all_59.bed 


