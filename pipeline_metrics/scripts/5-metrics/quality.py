import os
import sys
import pysam
import numpy as np
import math


pileup_file = sys.argv[1]
region_dir = sys.argv[2]
out_dir = sys.argv[3]

def main():
    pileup_tbx = pysam.TabixFile(pileup_file)
    
    # cycle through each gene in region file
    for region_file in os.listdir(region_dir):
        region_bed = open(region_dir+region_file, 'r')
        for line in region_bed:
            line = line.rstrip().split()
            try: line_chrm = int(line[0])
            except ValueError: line_chrm = line[0]
            line_start = int(line[1])
            line_end = int(line[2])
            
            for pileup_hit in pileup_tbx.fetch(line_chrm, line_start, line_end, parser=pysam.asBed()):
                pileup_chrm = pileup_hit[0]
                pileup_start = pileup_hit[1]
                pileup_end = pileup_hit[2]
                pileup_base_read = pileup_hit[6] 
                
                quality_score = calculate_quality(pileup_base_read)
                print(str(pileup_chrm) + '\t' + str(pileup_start) + '\t' + str(pileup_end) +
                            '\t' + str(quality_score) + '\t' + str(region_file.split('.')[0]))

                
def calculate_quality(pileup_base_read):
    row_quality = []
    for q in pileup_base_read:
        ascii_quality = float(ord(q) - ord('!'))
        phred_score = 10. ** (-ascii_quality/10.)
        quality_again = -10 * math.log10(phred_score)
        row_quality.append(ascii_quality)
    return(np.average(row_quality))
main()
