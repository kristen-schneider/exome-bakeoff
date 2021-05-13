import os
import sys
import pysam

pileup_file = sys.argv[1]
region_dir = sys.argv[2]
out_dir = sys.argv[3]

print('hello world')

def main():
    pileup_tbx = pysam.TabixFile(pileup_file)

    # cycle through each gene in region file
    for region_file in os.listdir(region_dir):
        region_bed = open(region_dir+region_file, 'r')
        print(region_dir+region_file)
    

#main()
