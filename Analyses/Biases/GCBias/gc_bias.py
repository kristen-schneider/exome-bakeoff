import pysam
import numpy as np
import pandas as pd
import pickle
from os import path, makedirs
import matplotlib.pyplot as plt
import seaborn as sns
import math
import time
import sys
import re


sns.set(style="whitegrid")


"""
Command line arguments
1. Reference Genome (.fasta)
2. Bam File (.bam)
3. Regions File (.bed)
"""


"""
import os
os.chdir('Analyses/Biases/GCBias')
"""


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


"""
Give a fasta sequence header formatted like:
>2 dna:chromosome chromosome:GRCh37:2:1:243199373:1
Return the chromosome, the beginning position and the end position
"""


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


"""
Given the path to the reference sequence
Break the reference into 100bp windows
Calc the GC % of each window
Return a pandas DataFrame each rowing representing a 100bp window and columns formatted as:
Chromosome  Start   End     GC%
"""


def load_hg37_gc_df(ref):
    # strip the path
    filename = re.sub('.*/', '', ref)
    # strip the file extension
    extension_less_name = re.sub('\\.\\w*', '', filename)
    cache_name = 'cache_' + extension_less_name + '.pickle'
    if not path.exists(cache_name):
        print('Generating Reference % GC info')
        chromosomes = []
        start = []
        end = []
        gc_percent = []
        with open(ref) as fp:
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
        pickle.dump(hg37gc_df, open(cache_name, 'wb'))
        return hg37gc_df
    else:
        print('Loading Reference Cache % GC info')
        return pickle.load(open(cache_name, 'rb'))


"""
Given a pandas DataFrame output like the output of load_hg37_gc_df()
Create a histogram with a density plot overlay of the number of 100 base windows found at each % GC bin
"""


def plot_hg37_gc_bias(df):
    data = df[df['GC%'] != 0]
    x = data['GC%']
    ax = sns.distplot(x)
    plt.savefig('hg37_gc_dist.png')


"""
Give an df formated like the output of load_hg37_gc_df
Return a dictionary with chromosomes as keys and a list containing 2 things as a value
    0. a df of just that chromosome indexed by start location
    1. a df of just that chromosome indexed by end location
"""


def create_dict_of_indexed_dfs(df):
    chromosomes = [x for x in list(set(df['Chromosome'])) if x[0] != 'G' and x[0] != 'M']
    dict_of_indexed_dfs = {}
    for c in chromosomes:
        temp_start = df[df['Chromosome'] == c]
        temp_end = df[df['Chromosome'] == c]
        temp_start.index = temp_start['Start']
        temp_end.index = temp_start['End']
        dict_of_indexed_dfs[c] = [temp_start, temp_end]
    return dict_of_indexed_dfs


"""
Given 
1. dictionary formted like the output of create_dict_of_indexed_dfs, 
2. the chromosome (string)
3. start locations of the read (int)
4. end locations of the read (int)
Return a list containing 1 or 2 GC% (int)
"""


def get_gc_content(indexed_dic,chrom,start,end):
    start_key = int(re.sub('\\w\\w$','01',str(start)))
    end_key = int(re.sub('\\w\\w$', '01', str(start))) + 100
    start_row = None
    end_row = None
    try:
        start_row = indexed_dic[chrom][0][start_key]
    except KeyError:
        start_row = indexed_dic[chrom][0][(start >= indexed_dic[chrom][0]['Start']) & (start < indexed_dic[chrom][0]['End'])]

    try:
        end_row = indexed_dic[chrom][1][end_key]
    except KeyError:
        end_row = indexed_dic[chrom][1][(end >= indexed_dic[chrom][1]['Start']) & (end < indexed_dic[chrom][1]['End'])]

    if int(end_row['Start']) != int(start_row['Start']):
        return [int(start_row['GC%'] * 100), int(end_row['GC%']* 100)]
    else:
        return [int(start_row['GC%'] * 100)]


"""
Given a reference genome file path, a bam file, a bed (or regions file) and the path for where results should be placed
Filter the bam to only use portions included in the bed file
Calculate normalized coverage across the reference genome % GC bins
Outputs a tab separated file to results_dir/GCBias/NormalizedCoverage/<name of bam file>.tsv
"""


def calc_gc_bias(ref, bam_file_path, bed, results_dir):
    # load the reference genome broken up into 100bp section annotated with gc %
    df = load_hg37_gc_df(ref)
    index_dict = create_dict_of_indexed_dfs(df)
    start = time.time()
    gc_bins_coverage = np.zeros(100)
    number_of_secondary_alignment_reads = 0

    bam = pysam.AlignmentFile(bam_file_path, "rb")
    for line in open(bed, 'r'):
        print(line)
        bed_line = line.split('\t')
        contig = bed_line[0]
        start = int(bed_line[1])
        end = int(bed_line[2])
        for f in bam.fetch(contig=contig, start=start, stop=end, until_eof=True):
            if f.is_secondary:
                number_of_secondary_alignment_reads += 1
                continue
            begin = f.pos
            width = len(f.seq)
            end = begin + width
            # print(str(begin) + ' ' + str(end))
            chr = f.rname
            if chr == 22:
                chr = 'X'
            elif chr == 23:
                chr = 'Y'
            else:
                chr = str(chr + 1)
            gc_list = get_gc_content(index_dict, chr, begin, end)
            for gc in gc_list:
                gc_bins_coverage[gc] += 1

    print('Duration: ' + str(time.time() - start))
    print('Number of secondary alignments excluded: ' + str(number_of_secondary_alignment_reads))

    # calculate normalized coverage (normalized to the mean)
    gc_mean = np.sum(gc_bins_coverage) / len(gc_bins_coverage)
    norm_gc = gc_bins_coverage / gc_mean

    # create df of the GC bias data
    gc_df = pd.DataFrame({'GC%': list(range(1, 101)), 'Coverage': gc_bins_coverage, 'Normalized_coverage': norm_gc})

    # form the output file path
    # strip the path
    filename = re.sub('.*/', '', bam_file_path)
    # strip the file extension2
    extension_less_name = re.sub('\\.\\w*', '', filename)
    if results_dir[-1] != '/':
        results_dir += '/'
    makedirs(results_dir)
    makedirs(results_dir + 'GCBias/')
    makedirs(results_dir + 'GCBias/NormalizedCoverage/')
    outfile = results_dir + 'GCBias/NormalizedCoverage/' + extension_less_name + '.tsv'
    gc_df.to_csv(outfile, sep='\t')


if __name__ == "__main__":
    ref = sys.argv[1]
    bam = sys.argv[2]
    bed = sys.argv[3]
    results_dir = sys.argv[4]
    calc_gc_bias(ref, bam, bed, results_dir)

9