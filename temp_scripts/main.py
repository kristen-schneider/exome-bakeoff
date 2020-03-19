import sys
import os
from get_metrics import gm_main
from plot_metrics import *

def main():
    pileup_path = sys.argv[1]
    region_path = sys.argv[2]
    final_metrics_path = '/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/chco-exome-analysis/final_metrics'    

    hm_Q_txt = open('heat_map_quality.txt', 'a')
    hm_Q_txt.truncate(0)
    hm_Q_txt.close()
    
    hm_SB_txt = open('heat_map_strand_bias.txt', 'a')
    hm_SB_txt.truncate(0)
    hm_SB_txt.close()
    
    num_samples = len([s for s in os.listdir(pileup_path) if s.endswith('bed.gz')])
    print (num_samples)
    for region_file in os.listdir(region_path):
        gm_main(pileup_path, region_path+region_file)
        pm_main(final_metrics_path, region_file, num_samples)
    
    plot_heat_map(num_samples)
main()  
