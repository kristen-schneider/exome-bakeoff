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
import matplotlib.patches as mpatches

sns.set_style("whitegrid", {'axes.grid': False})

"""
Given a DNA sequencing (assuming in all capital letters)
Return the % GC content
"""


def calc_percent_gc(seq):
    atgc_count = seq.count('C') + seq.count('T') + seq.count('A') + seq.count('G')
    if atgc_count == 0:
        return 0
    gc = (seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')) / atgc_count
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


def get_ref_region_df(ref, bed):
    filename = re.sub('.*/', '', ref)
    # strip the file extension
    extension_less_name = re.sub('\\.\\w*', '', filename)
    cache_name = 'cache_' + extension_less_name + '.pickle'
    if path.exists(cache_name) and False:
        print('Loading Reference and Bed Cache % GC info')
        return pickle.load(open(cache_name, 'rb'))
    print('Generating Reference and Bed Cache % GC info')
    bed_df = pd.read_csv(bed, sep='\t', header=None)
    bed_df.iloc[:, 0] = [str(x) for x in bed_df.iloc[:, 0]]
    pickle.dump(bed_df, open('bed_df.pickle', 'wb'))
    pickle.dump(ref, open('ref.pickle', 'wb'))
    bed_gc_bin_counts = np.zeros(101)
    ref_gc_bin_counts = np.zeros(101)
    individual_gene_counts_dict = {}
    # read in genome one chromosome / line at a time
    for name, seq in read_fasta(open(ref)):
        c, b, e = get_chrom_begin_end(name)
        # get the info from the bed that is part of the current chromosome
        sub_df = bed_df[bed_df.iloc[:, 0] == c]

        # break the reference into 100 bp section and bin them based on gc content
        for i in range(0, len(seq), 100):
            sub_string = seq[i: i + 100]
            dec_gc = calc_percent_gc(sub_string)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            gc = int((dec_gc + .001) * 100)
            ref_gc_bin_counts[gc] += 1

        # check if the bed needs any information from this chromosome
        if sub_df.shape[0] == 0:
            continue

        # get a strings of each portion of the bed and calc the gc content of that string
        for i in range(sub_df.shape[0]):
            start = sub_df.iloc[i, 1]
            end = sub_df.iloc[i, 2]
            gene = sub_df.iloc[i, -1]
            if gene not in individual_gene_counts_dict:
                individual_gene_counts_dict[gene] = np.zeros((101))
            string_start = b + start
            string_end = b + end
            # break the string into 100 bp sections
            for j in range(string_start, string_end, 100):
                sub_string = seq[j: j + 100]
                dec_gc = calc_percent_gc(sub_string)
                # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
                gc = int((dec_gc + .001) * 100)
                bed_gc_bin_counts[gc] += 1
                individual_gene_counts_dict[gene][gc] += 1

    df = pd.DataFrame(
        {'percent': list(range(0, 101)), 'ref_count': ref_gc_bin_counts, '59 genes count': bed_gc_bin_counts})
    for gene in individual_gene_counts_dict:
        df[gene] = individual_gene_counts_dict[gene]
    pickle.dump(df, open(cache_name, 'wb'))
    return df


# get the count for the BAMs
def calc_gc(bam, bed):
    percent_count = list(np.zeros((101)))
    # with that bed file, calc the % GC for all reads in the bam
    bam_obj = pysam.AlignmentFile(bam, "rb", require_index=False)
    for line in open(bed, 'r'):
        bed_line = line.split('\t')
        bed_contig = bed_line[0]
        bed_start = int(bed_line[1])
        bed_end = int(bed_line[2])
        for read in bam_obj.fetch(contig=bed_contig, start=bed_start, stop=bed_end, until_eof=True):
            gc_content = calc_percent_gc(read.seq)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            index = int((gc_content + .001) * 100)
            percent_count[index] += 1
    return percent_count


# get the count for the BAMs
def calc_gc_per_site(bam, bed):
    info = {'chr': [], 'start': [], 'end': [], 'avg_gc': [], 'gene': []}
    # with that bed file, calc the % GC for all reads in the bam
    bam_obj = pysam.AlignmentFile(bam, "rb", require_index=False)
    for line in open(bed, 'r'):
        bed_line = line.split('\t')
        bed_contig = bed_line[0]
        bed_start = int(bed_line[1])
        bed_end = int(bed_line[2])
        gene = bed_line[-1].strip()
        count = 0
        gc = 0
        for read in bam_obj.fetch(contig=bed_contig, start=bed_start, stop=bed_end, until_eof=True):
            gc_content = calc_percent_gc(read.seq)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            index = int((gc_content + .001) * 100)
            gc += index
            count += 1
        info['chr'].append(bed_contig)
        info['start'].append(bed_start)
        info['end'].append(bed_end)
        info['gene'].append(gene)
        # make sure there were some reads at this position
        if count != 0:
            info['avg_gc'].append(gc / count)
        else:
            info['avg_gc'].append('N\A')
    df = pd.DataFrame(info)
    return df


def get_bams_gc_df(bams, bed):
    sample_percents = {}
    for file in bams:
        if file[-4:] == '.bam':
            head, tail = os.path.split(file)
            tail = tail.replace('.bam', '')
            percents = calc_gc(file, bed)
            sample_percents[tail] = percents
    return pd.DataFrame(sample_percents)


def get_gc_bias_df(ref, bed, bams):
    bam_df = get_bams_gc_df(bams, bed)
    ref_reg_df = get_ref_region_df(ref, bed)
    pickle.dump(ref_reg_df, open('ref_reg_df.pickle', 'wb'))
    merged_df = pd.concat([ref_reg_df, bam_df], axis=1).fillna(0)
    return merged_df, ref_reg_df.shape[1], bam_df.shape[1]


def make_line_plot(df, ref_col, region_col, bam_indexes, sample_names, y_lab, fig_name, use_ref_and_reg=True):
    df = df[df['percent'] > 0]
    samp_name_set = list(set(sample_names))
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']
    widths = list(range(len(samp_name_set), 0, -1))
    samples = [list(df.columns)[x] for x in bam_indexes]
    # create dictionaries of the group and the color / width of their lines
    color_map = {samp_name_set[i]: colors[i] for i in range(len(samp_name_set))}
    width_map = {samp_name_set[i]: widths[i] for i in range(len(samp_name_set))}
    custom_lines = []

    fig, ax = plt.subplots()
    fig.set_size_inches(17.5, 12.5)
    samp_names_ordered = []
    for i in range(len(samples)):
        sample = samples[i]
        samp_name = sample_names[i]
        if samp_name not in samp_names_ordered:
            samp_names_ordered.append(samp_name)
            custom_lines.append(mpatches.Patch(color=color_map[samp_name]))
        ax.plot(df['percent'], df[sample], color=color_map[samp_name], linewidth=width_map[samp_name])

    if use_ref_and_reg:
        ax2 = ax.twinx()
        plotting_df = pd.DataFrame({'ref_col': list(df[ref_col]), 'reg_col': list(df[region_col])})
        plotting_df.plot.line(ax=ax2)
        ax2.set_xlabel('Percent GC')
        ax2.set_ylabel('Normalized Coverage')
        ax2.legend(loc='upper right', frameon=False)

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_ylabel(y_lab)

    # plt.legend(custom_lines, samples, frameon=False, loc='upper left')
    ax.legend(custom_lines, samp_names_ordered, frameon=False, loc='upper left')
    plt.xlabel('% GC')
    plt.savefig(fig_name)
    plt.clf()


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

    # build the df of all data
    df, ref_shape, bam_shape = get_gc_bias_df(ref, combine_bed_path, bams)
    df['ref_norm'] = df['ref_count'] / sum(df['ref_count'])
    df['region_norm'] = df['59 genes count'] / sum(df['59 genes count'])
    print(df.columns)
    quit()
    num_previously_added_cols = 2
    pickle.dump(df, open('samples_run_analyses.pickle', 'wb'))
    # split the sample named into their library prep and capture tech
    bam_indexes = list(range(ref_shape, ref_shape + bam_shape))
    sample_names = [list(df.columns)[x] for x in bam_indexes]
    tech_one = [re.sub('-.*-.*', '', x) for x in sample_names]
    tech_two = [re.sub('.*-(\\w+)-.*', '\\1', x) for x in sample_names]

    reg_df = df.iloc[:, 3:ref_shape]
    samp_df = df.iloc[:, ref_shape:ref_shape + bam_shape]

    make_heat_map_and_related_plots_for_individual_genes(reg_df, sample_names, combine_bed_path, bams, results_dir)
    make_heat_map_and_related_plots(reg_df, samp_df, results_dir, combine_bed_path)

    # normalize each of the samples
    norm_names = [[list(df.columns)[x], 'norm_' + list(df.columns)[x]] for x in bam_indexes]
    for names in norm_names:
        df[names[1]] = df[names[0]] / sum(df[names[0]])
    norm_bam_indexes = [x + len(bam_indexes) + num_previously_added_cols for x in bam_indexes]
    num_previously_added_cols += len(bam_indexes)

    # make joint probability of each column
    joint_names = [['norm_' + list(df.columns)[x], 'joint_' + list(df.columns)[x]] for x in bam_indexes]
    for names in joint_names:
        df[names[1]] = df[names[0]] * df['region_norm']
    joint_prob_bam_indexes = [x + len(bam_indexes) + num_previously_added_cols for x in bam_indexes]
    num_previously_added_cols += len(bam_indexes)

    # get expected coverage of each sample
    exp_names = [[list(df.columns)[x], 'exp_' + list(df.columns)[x]] for x in bam_indexes]
    for names in exp_names:
        df[names[1]] = df['region_norm'] * sum(df[names[0]])
    exp_bam_indexes = [x + len(bam_indexes) + num_previously_added_cols for x in bam_indexes]
    num_previously_added_cols += len(bam_indexes)

    # observed vs expected
    o_e_dif_names = [[list(df.columns)[x], 'o_e_dif_' + list(df.columns)[x]] for x in bam_indexes]
    for names in o_e_dif_names:
        df[names[1]] = df[names[0]] - df['exp_' + names[0]]
    ob_exp_bam_indexes = [x + len(bam_indexes) + num_previously_added_cols for x in bam_indexes]
    num_previously_added_cols += len(bam_indexes)

    # make the line plots
    make_line_plot(df, 'ref_norm', 'region_norm', bam_indexes, tech_one, 'Count',
                   results_dir + 'raw_count_library_prep.png')
    make_line_plot(df, 'ref_norm', 'region_norm', bam_indexes, tech_two, 'Count',
                   results_dir + 'raw_count_capture_tech.png')
    # normalized count
    make_line_plot(df, 'ref_norm', 'region_norm', norm_bam_indexes, tech_one, 'Normalized count',
                   results_dir + 'norm_count_library_prep.png')
    make_line_plot(df, 'ref_norm', 'region_norm', norm_bam_indexes, tech_two, 'Normalized count',
                   results_dir + 'norm_count_capture_tech.png')
    # joint probability
    make_line_plot(df, 'ref_norm', 'region_norm', joint_prob_bam_indexes, tech_one, 'Joint probability',
                   results_dir + 'joint_prob_library_prep.png')
    make_line_plot(df, 'ref_norm', 'region_norm', joint_prob_bam_indexes, tech_two, 'Joint probability',
                   results_dir + 'joint_prob_capture_tech.png')
    # observed expected difference
    make_line_plot(df, 'ref_norm', 'region_norm', ob_exp_bam_indexes, tech_one, 'Observed - Expected',
                   results_dir + 'obs_exp_library_prep.png', False)
    make_line_plot(df, 'ref_norm', 'region_norm', ob_exp_bam_indexes, tech_two, 'Observed - Expected',
                   results_dir + 'obs_exp_capture_tech.png', False)


"""
Given:
regions_df - (pandas DataFrame) columns are genes, rows are percent GC, cells are coverage at that % gc
samples_df - (pandas DataFrame) columns are samples, rows are percent GC, cells are coverage at that % gc
Computes an SSD (sum squared difference) of observed vs expected % gc coverage for each sample and each gene
Returns pandas DataFrame, list of column names, list of row names
"""


def get_ssd_heat_map_data(regions_df, samples_df):
    pickle.dump(regions_df, open('regions_df_get_ssd_heat_map_data.pickle', 'wb'))
    pickle.dump(samples_df, open('samples_df_get_ssd_heat_map_data.pickle', 'wb'))
    ssd_dict = {}
    tidy = {}
    tidy['sample'] = []
    tidy['gene'] = []
    tidy['ssd'] = []
    tidy['library_prep'] = []
    tidy['capture_tech'] = []
    for reg_col in list(regions_df.columns):
        reg_sub = regions_df[reg_col]
        reg_sub_norm = reg_sub / sum(reg_sub)
        gene_ssd = []
        for samp_col in list(samples_df.columns):
            expected = reg_sub_norm * sum(samples_df[samp_col])
            diff = expected - samples_df[samp_col]
            ssd = sum(diff * diff)
            gene_ssd.append(ssd)
            lib_prep = re.sub('-.*-.*', '', samp_col)
            capture_tech = re.sub('.*-(\\w+)-.*', '\\1', samp_col)
            tidy['sample'].append(samp_col)
            tidy['gene'].append(reg_col)
            tidy['ssd'].append(ssd)
            tidy['library_prep'].append(lib_prep)
            tidy['capture_tech'].append(capture_tech)
        ssd_dict[reg_col] = gene_ssd
    df = pd.DataFrame(ssd_dict)
    df.index = list(samples_df.columns)
    tidy_df = pd.DataFrame(tidy)

    return df, list(df.columns), list(df.index)


"""
Given:
df - (pandas DataFrame) columns are genes, rows are samples
tech_col - (string), name of column to be removed from dataframe and used as technology color labels
title - (string) title of plot
figname - (string) path to where the figure should be saved
colorscheme - (string) default 'RdBu' the matplotlib color scheme to be used
Creates hierarchically clustered heat map of the data and saves at the provided location
"""


def plot_clustered_heat_map(df, tech_col, title, figname, colorscheme='RdBu'):
    technology = df.pop(tech_col)

    lut = dict(zip(np.unique(technology), 'mkygbr'))
    row_colors = pd.Series(technology, name='library prep tech').map(lut)

    g = sns.clustermap(df,
                       figsize=(25, 30),
                       cmap=colorscheme,
                       row_colors=row_colors)

    g.fig.suptitle(title)

    plt.savefig(figname,
                dpi=150,
                figsize=(25, 30))


"""
Given:
df - (pandas DataFrame) columns are genes, rows are samples
lib_col - (string), name of column to be removed from dataframe and used as library_prep color labels
capture_col - (string), name of column to be removed from dataframe and used as capture_technology color labels
title - (string) title of plot
figname - (string) path to where the figure should be saved
colorscheme - (string) default 'RdBu' the matplotlib color scheme to be used
Creates hierarchically clustered heat map of the data and saves at the provided location
"""


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


"""
Given:
regions_df - (pandas DataFrame) columns are genes, rows are percent GC, cells are coverage at that % gc
samples_df - (pandas DataFrame) columns are samples, rows are percent GC, cells are coverage at that % gc
results_dir - (string) path to where figures should be saved
Computes an SSD (sum squared difference) of observed vs expected % gc coverage for each sample and each gene
Plots SSD using hierarchically save two figures, one for library prep, one for technology
"""


def make_heat_map_and_related_plots(reg_df, samp_df, results_dir, bed):
    ssd, col_names, row_names = get_ssd_heat_map_data(reg_df, samp_df)

    bed_df = pd.read_csv(bed, sep='\t', header=None)

    plot_ssd_vs_gene_length(ssd, bed_df, results_dir)

    ssd_both = ssd.copy()
    ssd_both['lib_prep_technology'] = [re.sub('-.*-.*', '', x) for x in list(ssd_both.index)]
    ssd_both['capture_technology'] = [re.sub('.*-(\\w+)-.*', '\\1', x) for x in list(ssd_both.index)]
    named_indexes = [list(ssd_both['lib_prep_technology'])[i] + ' ' + list(ssd_both['capture_technology'])[i] + ' ' * i
                     for i in range(ssd_both.shape[0])]
    ssd_both.index = named_indexes

    plot_clustered_heat_map_double_label(ssd_both, 'lib_prep_technology', 'capture_technology',
                                         'SSD of expected vs observed % GC coverage',
                                         results_dir + 'combine_ssd_heat_map.png', 'Blues_r')


def calc_gc_from_bed_df(bam, bed_df):
    percent_count = list(np.zeros((101)))
    # with that bed file, calc the % GC for all reads in the bam
    bam_obj = pysam.AlignmentFile(bam, "rb", require_index=False)
    for i in range(bed_df.shape[0]):
        bed_contig = bed_df.iloc[i, 0]
        bed_start = int(bed_df.iloc[i, 1])
        bed_end = int(bed_df.iloc[i, 2])
        for read in bam_obj.fetch(contig=bed_contig, start=bed_start, stop=bed_end, until_eof=True):
            gc_content = calc_percent_gc(read.seq)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            index = int((gc_content + .001) * 100)
            percent_count[index] += 1
    return percent_count


def make_heat_map_and_related_plots_for_individual_genes(reg_df, sample_names, bed_path, bams, results_dir):
    bed_df = pd.read_csv(bed_path, sep='\t', header=None)
    gene_dict = {}
    for gene in set(bed_df[4]):
        sub = bed_df[bed_df[4] == gene]
        gene_ssds = []
        gene_ref_coverages = reg_df[gene]
        gene_ref_total = sum(gene_ref_coverages)
        for bam in bams:
            observed = calc_gc_from_bed_df(bam, sub)
            total = sum(observed)
            differences = [observed[i] - (total * (gene_ref_coverages[i] / gene_ref_total)) for i in
                           range(len(observed))]
            ssd = sum([x * x for x in differences])
            gene_ssds.append(ssd)
        gene_dict[gene] = gene_ssds
    ssds = pd.DataFrame(gene_dict)
    ssds.index = sample_names
    ssd_both = ssds.copy()
    ssd_both['lib_prep_technology'] = [re.sub('-.*-.*', '', x) for x in list(ssd_both.index)]
    ssd_both['capture_technology'] = [re.sub('.*-(\\w+)-.*', '\\1', x) for x in list(ssd_both.index)]
    named_indexes = [list(ssd_both['lib_prep_technology'])[i] + ' ' + list(ssd_both['capture_technology'])[i] + ' ' * i
                     for i in range(ssd_both.shape[0])]
    ssd_both.index = named_indexes

    plot_clustered_heat_map_double_label(ssd_both, 'lib_prep_technology', 'capture_technology',
                                         'SSD of expected vs observed % GC coverage',
                                         results_dir + 'combine_ssd_heat_map_individual_genes.png', 'Blues_r')


def plot_ssd_vs_gene_length(ssd, bed_df, results_dir):
    # scatter plot
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    gene_lengths = {}
    ssd_means = {}
    for gene in ssd.columns:
        sub = bed_df[bed_df[4] == gene]
        # get the length, each end - start summed together
        length = sum(sub[2] - sub[1])
        gene_lengths[gene] = length
        ssd_means[gene] = sum(ssd[gene]) / ssd.shape[1]
        plt.scatter(([length] * ssd.shape[0]), ssd[gene])
    plt.ylabel('Observed vs expected GC content SSD')
    plt.xlabel('Log gene length')
    ax.set_xscale('log')
    plt.savefig(results_dir + 'gc_ssd_vs_gene_length.png')
    plt.clf()
    cor_coe = np.corrcoef(list(ssd_means.values()), list(gene_lengths.values()))[0, 1]
    with open(results_dir + 'gene_size_ssd_correlation.txt', 'w') as file:
        file.write('Pearson product-moment correlation coefficient: ' + str(cor_coe))
    # box plots
    # sort the genes by their length
    sorted_gene_lengths = {k: v for k, v in sorted(gene_lengths.items(), key=lambda item: item[1])}
    sorted_gene_names = list(sorted_gene_lengths.keys())
    ssd_sorted = ssd.reindex(sorted_gene_names, axis=1)
    fig, ax = plt.subplots()
    fig.set_size_inches(15, 10)
    ssd_sorted.plot(kind='box', ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
    plt.ylabel('Observed vs expected GC content SSD')
    plt.xlabel('Genes, sorted by ascending length')
    plt.savefig(results_dir + 'gc_ssd_vs_gene_length_boxplots.png')

def get_ref_bed_gc_content(ref,bed):
    bed_df = pd.read_csv(bed, sep='\t', header=None)
    bed_df.iloc[:, 0] = [str(x) for x in bed_df.iloc[:, 0]]
    out_data = {'chr':[],'start':[],'end':[],'gene':[],'gc':[],}
    for name, seq in read_fasta(open(ref)):
        c, b, e = get_chrom_begin_end(name)
        sub_df = bed_df[bed_df.iloc[:, 0] == c]
        for i in range(sub_df.shape[0]):
            start = sub_df.iloc[i, 1]
            end = sub_df.iloc[i, 2]
            gene = sub_df.iloc[i, -1]
            chrom = sub_df.iloc[i, 0]
            # get the reference sequence for this region of the genome
            gene_seq = seq[start-b:end-b]
            gc_content = calc_percent_gc(gene_seq)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            gc_int = int((gc_content + .001) * 100)
            out_data['chr'].append(chrom)
            out_data['start'].append(start)
            out_data['end'].append(end)
            out_data['gene'].append(gene)
            out_data['gc'].append(gc_int)
    return pd.DataFrame(out_data)



def subtract(items1, items2):
    items = []
    for i in range(len(items1)):
        try:
            items.append(float(items1[i]) - float(items2[i]))
        except:
            items.append('N\A')
    return items


"""
This tool calculates the % GC for each genome range in a regions bed file
To use this as a commandline tool it requires 3 parameters:
1. path to bam file
2. path to regions file in bed format
3. name of output file
"""

if __name__ == "__main__":
    bam = sys.argv[1]
    bed = sys.argv[2]
    out = sys.argv[3]
    ref = sys.argv[4]
    ref_region_cache_name = bed.split('/')[-1] + ref.split('/')[-1] + '_ref_region_temp.csv'
    if os.path.exists(ref_region_cache_name) and False:
        print('Loading Ref GC DF')
        ref_df = pd.read_csv(ref_region_cache_name)
    else:
        print('Making Ref GC DF')
        ref_df = get_ref_bed_gc_content(ref,bed)
        ref_df.to_csv(ref_region_cache_name)
    gc_df = calc_gc_per_site(bam, bed)
    print(list(gc_df['avg_gc']))
    print(type(list(gc_df['avg_gc'])[0]))
    gc_df['diff_from_ref'] = subtract(list(gc_df['avg_gc']), list(ref_df['gc']))
    gc_df.to_csv(out, sep='\t', index=False)
