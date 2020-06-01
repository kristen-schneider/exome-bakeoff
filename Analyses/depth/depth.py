import pysam
import pickle
from os import path, makedirs, listdir
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import re
import numpy as np
import sys

sns.set_style("whitegrid", {'axes.grid': False})


def plot_clustered_heat_map_double_label(df, lib_col, capture_col, title, figname, colorscheme='RdBu'):
    lib_prep_technology = df.pop(lib_col)
    capture_technology = df.pop(capture_col)

    lut = dict(zip(np.unique(lib_prep_technology), 'mkygbr'))
    lib_row_colors = pd.Series(lib_prep_technology, name='library prep tech').map(lut)

    lut = dict(zip(np.unique(capture_technology), 'mkygbr'))
    capture_row_colors = pd.Series(capture_technology, name='capture tech').map(lut)

    combine_colors = pd.DataFrame(lib_row_colors).join(pd.DataFrame(capture_row_colors))

    g = sns.clustermap(df,
                       figsize=(25, 30),
                       cmap=colorscheme,
                       row_colors=combine_colors)

    g.fig.suptitle(title)

    plt.savefig(figname,
                dpi=150,
                figsize=(25, 30))


def get_gene_coverage(bams, bed):
    print(bed)
    bed_df = pd.read_csv(bed, sep='\t', header=None)
    gene_coverage_dict = {}

    gene_lengths = {}
    for gene in set(bed_df[4]):
        sub = bed_df[bed_df[4] == gene]
        # get the length, each end - start summed together
        length = sum(sub[2] - sub[1])
        gene_lengths[gene] = length

    bam_index = 0
    for bam in bams:
        head, tail = os.path.split(bam)
        tail = tail.replace('.bam', '')
        print(tail)
        gene_coverage_dict[tail] = []
        print(str(bam_index) + '/' + str(len(bams)))
        bam_index += 1
        bam_obj = pysam.AlignmentFile(bam, "rb", require_index=False)
        line_count = 0
        gene_coverages = {gene: 0 for gene in set(bed_df[4])}
        for line in open(bed, 'r'):
            bed_line = line.strip().split('\t')
            bed_contig = str(bed_line[0])
            bed_start = int(bed_line[1])
            bed_end = int(bed_line[2])
            bed_gene = str(bed_line[4])
            count = 0
            for read in bam_obj.fetch(contig=bed_contig, start=bed_start, stop=bed_end, until_eof=True):
                count += 1
                gene_coverages[bed_gene] += 1
        gene_coverage_dict[tail] = list(gene_coverages.values())

    coverage_df = pd.DataFrame(gene_coverage_dict).T
    coverage_df.columns = list(gene_lengths.keys())
    norm_df = coverage_df.copy()
    for col in norm_df.columns:
        norm_df[col] = norm_df[col] / gene_lengths[col]

    return coverage_df, norm_df


def run_analyses(ref, bams, beds, results_dir):
    # make location to store intermediate files
    if results_dir[-1] != '/':
        results_dir += '/'
    intermediate_files_path = results_dir + 'Intermediate-Files/'
    try:
        makedirs(intermediate_files_path)
    except FileExistsError:
        pass

    # combind all the beds
    print('combining bed files')
    combine_bed_path = intermediate_files_path + 'combine_regions.bed'
    combind_bed = open(combine_bed_path, 'w')
    for bed in beds:
        for line in open(bed, 'r'):
            combind_bed.write(line)
    combind_bed.close()

    lib_preps = []
    capt_tech = []
    for bam in bams:
        head, tail = os.path.split(bam)
        tail = tail.replace('.bam', '')
        lib_preps.append(re.sub('-.*-.*', '', tail))
        capt_tech.append(re.sub('.*-(\\w+)-.*', '\\1', tail))

    df, norm_df = get_gene_coverage(bams, combine_bed_path)
    print(df)
    print(len(lib_preps))
    print(lib_preps)
    df['lib_prep_technology'] = lib_preps
    df['capture_technology'] = capt_tech
    norm_df['lib_prep_technology'] = lib_preps
    norm_df['capture_technology'] = capt_tech
    named_indexes = [list(norm_df['lib_prep_technology'])[i] + ' ' + list(norm_df['capture_technology'])[i] + ' ' * i
                     for i in range(norm_df.shape[0])]
    norm_df.index = named_indexes
    df.index = named_indexes
    plot_clustered_heat_map_double_label(df, 'lib_prep_technology', 'capture_technology',
                                         'Depth of coverage',
                                         results_dir + 'depth.png', 'Blues')
    plot_clustered_heat_map_double_label(norm_df, 'lib_prep_technology', 'capture_technology',
                                         'Normalized depth of coverage',
                                         results_dir + 'normalized_depth.png', 'Blues')
    return df, norm_df


d = None
n = None
if __name__ == "__main__":
    ref = sys.argv[1]
    bam_template = sys.argv[2]
    bed_template = sys.argv[3]
    res = sys.argv[4]

    # ref = '/Users/michael/TESTBAMs/human_g1k_v37.fasta'
    # bam_template = '/Users/michael/TESTBAMs/'
    # bed_template = '/Users/michael/TESTBAMs/'
    # res = '/Users/michael/BakeOff/Results/GCBias/'

    # make there paths are not missing / on the end
    if bam_template[-1] != '/':
        bam_template += '/'
    if bed_template[-1] != '/':
        bed_template += '/'

    # make list of all bams
    bams = []
    for file in listdir(bam_template):
        if file[-4:] == '.bam':
            bams.append(bam_template + file)

    # make list of all beds
    beds = []
    for file in listdir(bed_template):
        if file[-4:] == '.bed':
            beds.append(bed_template + file)

    # run the analysis
    d, n = run_analyses(ref, bams, beds, res)
    pickle.dump(d, open('depth.pickle', 'wb'))
