import pysam
import numpy as np
import pandas as pd
import pickle
from os import path, makedirs
import seaborn as sns
import sys
import re
import multiprocessing as mp
import os
import pysamstats

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
    atgc_count = seq.count('C') + seq.count('T') + seq.count('A') + seq.count('G')
    gc = (seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')) / atgc_count
    at = (seq.count('A') + seq.count('T') + seq.count('a') + seq.count('t')) / atgc_count
    if gc + at != 1:
        print('FUCK! AT% ' + str(at) + ' GC% ' + str(gc) + ' length ' + str(len(seq))) + ' atgc len' + str(atgc_count)
        print(seq)
    return gc


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


def calc_gc_old(bam, bed):
    percent_count = list(np.zeros((101)))
    # with that bed file, calc the % GC for all reads in the bam
    bam_obj = pysam.AlignmentFile(bam, "rb", require_index=False)
    lengths = []
    for line in open(bed, 'r'):
        print(line)
        bed_line = line.split('\t')
        bed_contig = bed_line[0]
        bed_start = int(bed_line[1])
        bed_end = int(bed_line[2])
        for read in bam_obj.fetch(contig=bed_contig, start=bed_start, stop=bed_end, until_eof=True):
            gc_content = calc_percent_gc(read.seq)
            # lengths.append(len(read.seq))
            index = int(gc_content * 100)
            if index == 56:
                lengths.append(len(read.seq))
            percent_count[index] += 1

    head, tail = os.path.split(bam)

    df = pd.DataFrame({'percent': list(range(0, len(percent_count))), 'count': percent_count})
    df.to_csv('new-results/' + tail + '.csv')
    return df, lengths


def count_windows_at_percent(params):
    percent = params[0]
    print('starting percent ' + str(percent))
    df = params[1]
    per = percent / 100
    per_upper = per + .01
    df = df[df['GC%'] >= per]
    df = df[df['GC%'] < per_upper]
    print('finished percent ' + str(percent))
    return [percent, df.shape[0]]


def get_windows_count(ref):
    # produce the number of reference windows at each percentage
    windows = None
    windows_path = 'reference_windows_counts.pickle'
    if os.path.exists(windows_path) and False:
        print('Loading Windows Count')
        windows = pickle.load(open(windows_path, 'rb'))
    else:
        print('Generating Windows Count')
        ref_df = load_hg37_gc_df(ref)
        pool = mp.Pool(processes=12)
        reads_res = pool.map(count_windows_at_percent, [[i, ref_df] for i in range(0, 101)])
        window_percent = []
        window_count = []
        for i in range(len(reads_res)):
            window_percent.append(reads_res[i][0])
            window_count.append(reads_res[i][1])
        windows = pd.DataFrame({'percent': window_percent, 'window_count': window_count})
        pickle.dump(windows, open(windows_path, 'wb'))
    return windows


def calc_gc(ref, bam, bed, res):
    try:
        makedirs(res)
    except FileExistsError:
        pass

    gcs = []
    for line in open(bed, 'r'):
        # print(line)
        bed_line = line.split('\t')
        bed_contig = bed_line[0]
        bed_start = int(bed_line[1])
        bed_end = int(bed_line[2])
        for rec in pysamstats.stat_coverage_gc(bam, ref, chrom=bed_contig, start=bed_start, end=bed_end):
            gcs.append((rec['gc']))

    counts = np.zeros((101))
    for i in gcs:
        counts[i] += 1

    ref_windows = get_windows_count(start_ref)
    gc_df = pd.DataFrame({'percent': list(range(0, 101)), 'count': counts, 'ref_windows': ref_windows['window_count']})
    head, tail = os.path.split(bam)
    tail = tail.replace('.bam','')

    if res[-1] != '/':
        res += '/'
    gc_df.to_csv(res + tail + '.csv')

    return gc_df


if __name__ == "__main__":
    start_ref = sys.argv[1]
    start_bam = sys.argv[2]
    start_bed = sys.argv[3]
    start_res = sys.argv[4]
    calc_gc(start_ref, start_bam, start_bed, start_res)
