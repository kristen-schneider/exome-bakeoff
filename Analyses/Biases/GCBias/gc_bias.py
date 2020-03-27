import pysam
import numpy as np
import pandas as pd
import pickle
from os import path
import matplotlib.pyplot as plt
import seaborn as sns
import math
import time

sns.set(style="whitegrid")


"""
import os
os.chdir('Analyses/Biases/GCBias')
"""

# htslib
# set up a conda env with a specific version of pysam
# create a requirement.txt and list the version of everything we have used

# read in reference genome
# break reference into 100bp windows
# calc % GC for each window
# pickle this information

# df or matrix (np array?)
# X = 23, Y = 24, MT = 25, scaffold = 26
# chromosome    start   end     %GC


"""
Given a DNA sequencing (assuming in all capital letters)
Return the % GC content
"""


def calc_percent_gc(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)


"""
Taken from: https://stackoverflow.com/questions/29805642/learning-to-parse-a-fasta-file-with-python
Give an file object
return the sequence name and a single string of the bases
"""


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def get_chrom_begin_end(name):
    split_name = name.split(' ')
    chrom = split_name[0].replace('>', '')
    # control for mitochrondrial DNA what is formated different
    if chrom[0] == 'M':
        return chrom, 1, -1

    split_location = split_name[2].split(':')
    begin = int(split_location[-3])
    fin = int(split_location[-2])
    return chrom, begin, fin


def load_hg37_gc_df():
    pickle_path = 'exome_bakeoff_hg37_gc_pd.pickle'
    if not path.exists(pickle_path):
        print('Generating HG 37 % GC info')
        chromosomes = []
        start = []
        end = []
        gc_percent = []
        with open("/Users/michael/TESTBAMs/human_g1k_v37.fasta") as fp:
            for name, seq in read_fasta(fp):
                print(name)
                c, b, e = get_chrom_begin_end(name)
                print(str(c) + ' ' + str(b) + ' ' + str(e))
                gc_percent = gc_percent + [calc_percent_gc(seq[i:i + 100]) for i in range(0, len(seq), 100)]
                new_start = np.array(list(range(0, len(seq), 100)))
                new_start = new_start + b
                new_end = new_start + 100
                new_chrom = [c] * len(new_start)
                new_start.tolist()
                new_end.tolist()
                start = start + new_start.tolist()
                end = end + new_end.tolist()
                chromosomes = chromosomes + new_chrom
        hg37gc_df = pd.DataFrame({'Chromosome': chromosomes, 'Start': start, 'End': end, 'GC%': gc_percent})
        pickle.dump(hg37gc_df, open(pickle_path, 'wb'))
        return hg37gc_df
    else:
        print('Loading HG 37 % GC info')
        return pickle.load(open(pickle_path, 'rb'))


def plot_hg37_gc_bias(df):
    data = df[df['GC%'] != 0]
    x = data['GC%']
    ax = sns.distplot(x)
    plt.savefig('hg37_gc_dist.png')


df = load_hg37_gc_df()
# plot_hg37_gc_bias(df)


bam = pysam.AlignmentFile("/Users/michael/TESTBAMs/small-Twist-Roche-0720ME25_S5_L001.bqsr.bam", "rb")

# get the coverage of each window

sec = 0
pri = 0
thing = None

# filter reads down to only the primary alignment
thing = None
rnames = []
gc_bins_coverage = np.zeros(100)
chr_sub_set = None
current_chr = None
start = time.time()
error_count = 0
current = 0
current_count = 0
for r in bam:
    if r.is_secondary:
        # print('Secondary')
        sec += 1
        continue
    # this is for a speed up, so boolean opperations are not performed across all chromosomes, only the one in question
    # BAM files are sorted so this will only resolve true 23 times
    if current_chr is not r.rname:
        chr_sub_set = df[df['Chromosome'] == str(r.rname + 1)]
        current_chr = r.rname
        print('----' + str(r.rname) + '----')
        print(time.time() - start)
    # this is another speed up, every 10 thousand toss out the ones never to be needed again
    # if current == 10000:
    #     print('\t' + str(begin))
    #     chr_sub_set = chr_sub_set[chr_sub_set['Start'] > r.pos]
    #     current = 0
    # else:
    #     current += 1
    pri += 1
    begin = r.pos
    width = len(r.seq)
    # which bin the read beginning belongs in
    start_df = chr_sub_set[(begin >= chr_sub_set['Start']) & (begin < chr_sub_set['End'])]
    end_df = chr_sub_set[(begin+width >= chr_sub_set['Start']) & (begin+width < chr_sub_set['End'])]
    try:
        if int(start_df['Start']) == int(end_df['Start']):
            # the start and end fall in the same window of the genome
            start_percent = math.floor(start_df['GC%']*100)
            gc_bins_coverage[start_percent] += 1
        else:
            # the start and end are in two different regions
            start_percent = math.floor(start_df['GC%'] * 100)
            gc_bins_coverage[start_percent] += 1
            end_percent = math.floor(end_df['GC%'] * 100)
            gc_bins_coverage[end_percent] += 1
    except TypeError:
        error_count += 1
        print('----------------------------------------------------')
        print('Type Error, probably in the if statement casting to an int')
        print(start_df['Start'])
        print('---')
        print(end_df['Start'])
        print('----------------------------------------------------')


print('Error Count' + str(error_count))
print(time.time() - start)

# calcualte normalized coverage (normalized to the mean)
gc_mean = np.sum(gc_bins_coverage) / len(gc_bins_coverage)
norm_gc = gc_bins_coverage / gc_mean

gc_df = pd.DataFrame({'GC%':list(range(1,101)),'Coverage':gc_bins_coverage,'Normalized_coverage':norm_gc})

gc_df.to_csv('example.tsv', sep='\t')

print(sec)
print(pri)
