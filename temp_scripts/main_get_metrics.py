import sys
import os
from get_metrics import gm_main

def main():
    pileup_path = sys.argv[1]
    region_path = sys.argv[2]
    final_metrics_path = '/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/chco-exome-analysis/final_metrics'

    for region_file in os.listdir(region_path):
        gm_main(pileup_path, region_path+region_file)
main()
