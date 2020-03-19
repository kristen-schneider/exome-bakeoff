# imports
import sys
import os
import numpy as np
from scipy.stats import binom
import pysam
import csv
import matplotlib.pyplot as plt
from genome_to_exome_space import geno_to_exo_main

###########################
###########################
def gm_main(path, region_file):
    simple_file_name = ''
    simple_file_name = os.path.basename(region_file)
    simple_file_name = os.path.splitext(simple_file_name)[0]

    
    regions = get_chrm_start_end(region_file)
    #print(regions[0]
    # for each file in path
    for bed_file in os.listdir(path):
        
        if bed_file.endswith('bed.gz'):
            # file outputs
            # intermediate quality
            intermediates_path = '/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/chco-exome-analysis/intermediates'
            intermediate_Q_name = bed_file+simple_file_name+'quality_intermediate.txt'
            intermediate_Q_txt = open(os.path.join(intermediates_path, intermediate_Q_name), 'a')
            intermediate_Q_txt.truncate(0)
            intermediate_SB_name = bed_file+simple_file_name+'strand_bias_intermediate.txt'
            intermediate_SB_txt = open(os.path.join(intermediates_path, intermediate_SB_name), 'a')
            intermediate_SB_txt.truncate(0)

            # final quality
            finals_path = '/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/chco-exome-analysis/final_metrics'
            final_Q_name = bed_file+simple_file_name+'quality_final.txt'
            final_Q_txt = open(os.path.join(finals_path, final_Q_name), 'a')
            final_Q_txt.truncate(0)
            final_SB_name = bed_file+simple_file_name+'strand_bias_final.txt'
            final_SB_txt = open(os.path.join(finals_path, final_SB_name), 'a')       
            final_SB_txt.truncate(0)         


            for r in range(len(regions)):
                chrm = regions[r][0]
                start = int(regions[r][1])
                end = int(regions[r][2])
                tbx = pysam.TabixFile(path + '/' + bed_file)    
        
                # list to hold mpileup quality counts
                quality = get_quality(tbx, chrm, start, end)
                for q in quality: print(q[0], '\t', q[1], '\t', q[2], '\t',  np.average(q[3]), file=intermediate_Q_txt)
    
                reads = get_reads(tbx, chrm, start, end)
                counts = get_counts(reads)
                strand_bias = get_strandbias(counts)
                if strand_bias == -1: pass
                else:
                    for sb in strand_bias: print(sb[0][0], '\t', sb[0][1], '\t', sb[0][2], '\t', sb[1], file=intermediate_SB_txt)
            
            intermediate_Q_txt.close()
            intermediate_SB_txt.close()  
            geno_to_exo_main(region_file, os.path.join(intermediates_path, intermediate_Q_name), os.path.join(finals_path, final_Q_name))
            geno_to_exo_main(region_file, os.path.join(intermediates_path, intermediate_SB_name), os.path.join(finals_path, final_SB_name))
        
        else: continue
        #intermediate_Q_txt.close()
        #intermediate_SB_txt.close()

###########################
###########################
def get_chrm_start_end(region_file):
    region = []
    f = open(region_file, 'r')
    for line in f:
        current_region = []
        
        line = line.rstrip().split()
        chrm = line[0]
        start = line[1]
        end = line[2]
        current_region = [chrm, start, end]
        region.append(current_region)
    return region
    


###########################
# Input:
# 1. tbx mpileup file
# 2. chromosome
# 3. 1-start
# 4. open-end 
# Output: 
# [start, end, q_values (as given my mpileup notation)]
###########################
def get_quality(tbx, chrm, start, end):
    
    quality = []
    
    # pull out the quality scores
    for row in tbx.fetch(chrm, start, end, parser=pysam.asBed()):
        c = row[0]
        s = row[1]
        e = row[2]
        q = row[7]
        
        # [c, s, e, q_values]
        row_info = []
        row_info.extend((c, s, e))
        # to hold quality scores converted from ASCII from 0 to 93
        q_values = [] 
        for qv in q:
            q_values.append(ord(qv)-ord('!'))

        row_info.append(q_values)
        quality.append(row_info)
    
    return quality

###########################
###########################
def get_reads(tbx, chrm, start, end):
    
    reads = []
    
    # pull out the quality scores
    for row in tbx.fetch(chrm, start, end, parser=pysam.asBed()):
        c = row[0]
        s = row[1]
        e = row[2]
        r = row[5]
        cse = [c, s, e]
        curr_row = [cse, r]
        reads.append(curr_row)
    #for r in reads: print (r)
    return reads
            
    
###########################
###########################
def get_counts(reads):
    counts = [] 
    for r in reads:
        row_counts = []

        pos_match_count = 0
        pos_mismatch_count = 0
        neg_match_count = 0
        neg_mismatch_count = 0
       
        read_string = r[1]
        #print(read_string)
        for i in read_string:
           # print (i)
            if i == '.': pos_match_count += 1
            elif (i == 'A' or i == 'C' or i == 'T' or i == 'G'): pos_mismatch_count += 1
            elif i == ',': neg_match_count += 1
            elif (i == 'a' or i == 'c' or i == 't' or i == 'g'): neg_mismatch_count += 1
            else: continue
        row_counts = [pos_match_count, pos_mismatch_count, neg_match_count, neg_mismatch_count]
        full_row = [r[0], row_counts]
        #print(full_row)
        counts.append(full_row)
    return counts


###########################
###########################
def get_strandbias(counts):
    strand_bias = []
    
    for c_i in counts:
        curr_position = c_i[0]
        curr_counts = c_i[1]
        #print (c_i)
        a = float(curr_counts[0])
        b = float(curr_counts[1])
        c = float(curr_counts[2])
        d = float(curr_counts[3])

        try: sb = (c+d)/(a+b+c+d)
        except ZeroDivisionError:
            return -1
        curr_row = [curr_position, sb]
        strand_bias.append(curr_row)

    return strand_bias


#main ()
