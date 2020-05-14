import pysam
import pickle
from os import path, makedirs, listdir
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    if path.exists(cache_name):
        print('Loading Reference and Bed Cache % GC info')
        return pickle.load(open(cache_name, 'rb'))
    print('Generating Reference and Bed Cache % GC info')
    bed_df = pd.read_csv(bed, sep='\t')
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
    merged_df = pd.concat([ref_reg_df, bam_df], axis=1).fillna(0)
    return merged_df, ref_reg_df.shape[1], bam_df.shape[1]


def make_line_plot(df, ref_col, region_col, bam_indexes, sample_names, y_lab, fig_name, use_ref_and_reg=True):
    print(bam_indexes)
    df = df[df['percent'] > 0]
    samp_name_set = list(set(sample_names))
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']
    widths = list(range(len(samp_name_set), 0, -1))
    print(widths)
    print(df.shape)
    samples = [list(df.columns)[x] for x in bam_indexes]
    print(samples)
    # create dictionaries of the group and the color / width of their lines
    color_map = {samp_name_set[i]: colors[i] for i in range(len(samp_name_set))}
    width_map = {samp_name_set[i]: widths[i] for i in range(len(samp_name_set))}
    custom_lines = []

    fig, ax = plt.subplots()
    fig.set_size_inches(17.5, 12.5)
    for i in range(len(samples)):
        sample = samples[i]
        samp_name = sample_names[i]
        # custom_lines.append(Line2D([0], [0], color=color_map[sample], lw=4, label=sample + 'banene'))
        print(sample)
        print(df[sample])
        ax.plot(df['percent'], df[sample], color=color_map[samp_name], linewidth=width_map[samp_name])

    for samp_name in list(set(sample_names)):
        custom_lines.append(mpatches.Patch(color=color_map[samp_name]))
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
    ax.legend(custom_lines, sample_names, frameon=False, loc='upper left')
    plt.xlabel('% GC')
    plt.savefig(fig_name)
    plt.clf()


# make a heat map plotting function
def plot_heat_map(df, col_indexes, figname):
    columns = [list(df.columns)[x] for x in col_indexes]
    plot_data = df[df.columns.intersection(columns)]
    fig, ax = plt.subplots()
    fig.set_size_inches(35, 25)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns)
    ax.set_ylabel('% GC')
    plt.imshow(plot_data)
    plt.colorbar(cmap='cold')
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.savefig(figname)
    plt.clf()


def make_ssd_table(dd, figname):
    ssd_df = pd.DataFrame(dd).T
    ssd_df.columns = ['SSD', 'Library Prep', 'Capture Tech']
    # sort ascending
    ssd_df = ssd_df.sort_values('SSD')

    cells_ar = ssd_df.to_numpy()
    lib_preps = list(set(ssd_df['Library Prep']))
    capt_tech = list(set(ssd_df['Capture Tech']))
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']

    color_lib_map = {lib_preps[x]: colors[x] for x in range(len(lib_preps))}
    color_capt_map = {capt_tech[x]: colors[x] for x in range(len(capt_tech))}
    lib_colors = [color_lib_map[x] for x in list(ssd_df['Library Prep'])]
    tech_colors = [color_capt_map[x] for x in list(ssd_df['Capture Tech'])]
    background_df = ssd_df.copy()
    background_df['SSD'] = 'w'
    background_df['Library Prep'] = lib_colors
    background_df['Capture Tech'] = tech_colors

    fig, ax = plt.subplots()
    fig.set_size_inches(8.5, 15)
    ax.table(cellText=cells_ar, cellColours=background_df.to_numpy(), colLabels=['SSD', 'Library Prep', 'Capture Tech'],
             loc='center')
    ax.axis('tight')
    ax.axis('off')
    plt.savefig(figname)


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
    print(df)
    df['ref_norm'] = df['ref_count'] / sum(df['ref_count'])
    df['region_norm'] = df['59 genes count'] / sum(df['59 genes count'])
    num_previously_added_cols = 2

    # split the sample named into their library prep and capture tech
    bam_indexes = list(range(ref_shape, ref_shape + bam_shape))
    sample_names = [list(df.columns)[x] for x in bam_indexes]
    tech_one = [re.sub('-.*-.*', '', x) for x in sample_names]
    tech_two = [re.sub('.*-(\\w+)-.*', '\\1', x) for x in sample_names]
    pickle.dump(df,open('df.pickle','wb'))
    pickle.dump(tech_one, open('tech_one.pickle', 'wb'))
    pickle.dump(tech_two, open('tech_two.pickle', 'wb'))
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

    # SSD
    samp_ssd_dict = {}
    index = 0
    for x in bam_indexes:
        col_name = list(df.columns)[x]
        samp_ssd_dict[col_name] = []
        samp_ssd_dict[col_name].append(sum(df[col_name] * df[col_name]))
        samp_ssd_dict[col_name].append(tech_one[index])
        samp_ssd_dict[col_name].append(tech_two[index])
        index += 1

    # make table
    make_ssd_table(samp_ssd_dict, results_dir + 'ssd_table.png')
    # plot heat map
    plot_heat_map(df, norm_bam_indexes, results_dir + 'norm_heat_map.png')
    plot_heat_map(df, list(range(3, 61)), results_dir + '59_genes_heat_map.png')
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
    return df, norm_bam_indexes, samp_ssd_dict


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
    df, norm_indexes, samp_ssd_dict = run_analyses(ref, bams, beds, res)
    pickle.dump(samp_ssd_dict, open('samp_ssd_dict.pickle', 'wb'))
