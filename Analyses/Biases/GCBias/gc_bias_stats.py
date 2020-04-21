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

sns.set_style("whitegrid", {'axes.grid' : False})


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
Given
bam - path to a bam file
bed - path to a bed file
res - path of the results directory
Calculate the % GC of reads in the bam file that occur within the regions specified in the bed file
return a pandas DataFrame with columns ('percent','count')
"""


def calc_gc(bam, bed, res):
    percent_count = list(np.zeros((101)))
    # with that bed file, calc the % GC for all reads in the bam
    bam_obj = pysam.AlignmentFile(bam, "rb", require_index=False)
    for line in open(bed, 'r'):
        # print(line)
        bed_line = line.split('\t')
        bed_contig = bed_line[0]
        bed_start = int(bed_line[1])
        bed_end = int(bed_line[2])
        for read in bam_obj.fetch(contig=bed_contig, start=bed_start, stop=bed_end, until_eof=True):
            gc_content = calc_percent_gc(read.seq)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            index = int((gc_content + .001) * 100)
            percent_count[index] += 1

    head, tail = os.path.split(bam)
    df = pd.DataFrame({'percent': list(range(0, len(percent_count))), 'count': percent_count, })

    if res[-1] != '/':
        res += '/'
    tail = tail.replace('.bam', '')
    df.to_csv(res + tail + '.csv')
    return df


"""
Given
df - pandas DataFrame with columns ('percent','SAMPLE','tech_one','tech_two','ref_bin_normed_count') - not ordered
group_col -  the column to be used for grouping the data (tech_one or tech_two)
y_col - column to be used in the y axis
y_lab - label to be used on the y axis
fig_name - name of the figure to be produced
Creates a line plot based on the give parameters and saves it
"""


def plot_by_tech(df, group_col, y_col, y_lab, fig_name,scale=1):
    # rgba based color schemes to make colors semi transparent
    # pale = [(.07,.46,.2,.5), (.2,.14,.53,.5),(.86,.8,.46,.5),(.8,.4,.46,.5),(.53,.01,.33,.5),(.54,.8,.93,.5)]
    # opacity = .8
    # colors = [(230/255,159/255,3/255,opacity),(86/255,180/255,233/255,opacity),(5/255,158/255,115/255,opacity),(3/255,114/255,178/255,opacity),(213/255,94/255,0,opacity),(204/255,121/255,167/255,opacity)]
    colors = ['r', 'g', 'b', 'orange', 'purple', 'y']
    widths = [6, 5, 4, 3, 2, 1]
    samples = list(set(df[group_col]))
    # create dictionaries of the group and the color / width of their lines
    color_map = {samples[i]: colors[i] for i in range(len(samples))}
    width_map = {samples[i]: widths[i] for i in range(len(samples))}
    custom_lines = []

    ax = plt.subplot(111)
    for sample in samples:
        custom_lines.append(Line2D([0], [0], color=color_map[sample], lw=4))
        sub = df[df[group_col] == sample]
        for s in list(set(sub['SAMPLE'])):
            subsub = sub[sub['SAMPLE'] == s]
            subsub[y_col] = [x / scale for x in subsub[y_col]]
            ax.plot(subsub['percent'], subsub[y_col], color=color_map[sample],
                    linewidth=width_map[sample])

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.legend(custom_lines, samples, frameon=False, loc='upper left')
    plt.xlabel('% GC')
    plt.ylabel(y_lab)

    plt.savefig(fig_name, dpi=300)
    plt.clf()


"""
Given:
csv_dir - path to where the output csv from gc_bias_stats.py are stored
results_dir - path to directory the figures should be saved
"""


def make_plots(csv_dir, results_dir, ref_bin_counts):
    # add a back slash to the of the dirs in case there is no one already
    if results_dir[-1] != '/':
        results_dir += '/'
    if csv_dir[-1] != '/':
        csv_dir += '/'

    # make the results directory if it does not exist already
    try:
        makedirs(results_dir)
    except FileExistsError:
        pass

    data = None
    first = True
    first_df = None

    # read in all the gc bias stats .csv files and combine into one DataFrame
    first = True
    for f in os.listdir(csv_dir):
        if f[-3:] != 'csv':
            continue
        df = pd.read_csv(csv_dir + f)
        df['SAMPLE'] = f.replace('.csv', '')
        if first:
            data = df
            first = False
            first_df = df
        else:
            data = data.append(df)

    # create a column with the number of reads found in the reference genome at that region
    data['ref_windows'] = [ref_bin_counts[x] for x in list(data['percent'])]

    # create alternate columns to plot with
    data['ref_bin_normed_count'] = data['count'] / data['ref_windows']
    data['tech_one'] = [re.sub('-.*-.*', '', x) for x in list(data['SAMPLE'])]
    data['tech_two'] = [re.sub('.*-(\\w+)-.*', '\\1', x) for x in list(data['SAMPLE'])]

    plot_by_tech(data, 'tech_one', 'ref_bin_normed_count', 'Normalized Read Count',
                 results_dir + 'gc_bias_library_prep.png')
    plot_by_tech(data, 'tech_two', 'ref_bin_normed_count', 'Normalized Read Count',
                 results_dir + 'gc_bias_capture_technology.png')

    plot_by_tech(data, 'tech_one', 'count', 'Read count (thousands)', results_dir + 'raw_coverage_gc_bias_library_prep.png',1000)
    plot_by_tech(data, 'tech_two', 'count', 'Read count (thousands)', results_dir + 'raw_coverage_gc_bias_capture_technology.png',1000)


"""
Give the path to a .bed file and reference genome (.fasta)
Plot the coverage based on % GC of the whole reference genome and sections of the genome specified in the .bed
Returns a numpy array of length 101 (percents 0-100) representing the % GC coverage of the reference genome 
"""


def plot_59_genes_reference_gc_bias(bed, ref, res):
    bed_df = pd.read_csv(bed, sep='\t')
    bed_gc_bin_counts = np.zeros(101)
    ref_gc_bin_counts = np.zeros(101)
    # read in genome one chromosome / line at a time
    for name, seq in read_fasta(open(ref)):
        c, b, e = get_chrom_begin_end(name)
        # get the info from the bed that is part of the current chromosome
        sub_df = bed_df[bed_df.iloc[:, 0] == c]

        # check if the bed needs any information from this chromosome
        if sub_df.shape[0] == 0:
            continue

        # break the reference into 100 bp section and bin them based on gc content
        for i in range(0, len(seq), 100):
            sub_string = seq[i: i + 100]
            dec_gc = calc_percent_gc(sub_string)
            # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
            gc = int((dec_gc + .001) * 100)
            ref_gc_bin_counts[gc] += 1

        # get a strings of each portion of the bed and calc the gc content of that string
        for i in range(sub_df.shape[0]):
            start = sub_df.iloc[i, 1]
            end = sub_df.iloc[i, 2]
            string_start = b + start
            string_end = b + end
            # break the string into 100 bp sections
            for j in range(string_start, string_end, 100):
                sub_string = seq[j: j + 100]
                dec_gc = calc_percent_gc(sub_string)
                # + .001 is to account for floating point rounding errors in python 0.58 is really 57.99999999
                gc = int((dec_gc + .001) * 100)
                bed_gc_bin_counts[gc] += 1
    og_ref_gc_bin_counts = ref_gc_bin_counts.copy()
    ref_gc_bin_counts = [x / 1000 for x in ref_gc_bin_counts]

    df = pd.DataFrame(
        {'Percent': list(range(0, 101)), 'reference count': ref_gc_bin_counts, '59 genes count': bed_gc_bin_counts})
    df = df[df['Percent'] > 0]
    ax = df.plot(x="Percent", y="reference count", legend=False)
    ax2 = ax.twinx()
    df.plot(x="Percent", y="59 genes count", ax=ax2, legend=False, color="r")
    ax.figure.legend()
    ax.grid(False)
    ax2.grid(False)
    ax.set_ylabel('Reference genome count (thousands)')
    ax2.set_ylabel('59 Genes count')
    ax.set_xlabel('% GC')
    if res[-1] != '/':
        res += '/'
    plt.savefig(res + 'gc_bias_reference_and_59_genes.png')
    plt.clf()
    return og_ref_gc_bin_counts


"""
Given
ref - path to reference genome
bams - list of paths to bam files
beds - list of paths to bed files
results_dir - location results should be stored
Generates summary statistics for each of the bam files
Produced plots summarizing GC bias in the bams, reference genome and regions of interest
"""


def run_gc_bias_analysis(ref, bams, beds, results_dir):
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

    # create intermediate location for gc bias summary stats
    print('creating intermediate location')
    summary_stats_path = intermediate_files_path + 'Summary-Stats-CSVs/'
    try:
        makedirs(summary_stats_path)
    except FileExistsError:
        pass

    # generate gc bias stats for each bam
    print('doing GC stats for bams')
    for bam in bams:
        print('\t' + bam)
        calc_gc(bam, bed, summary_stats_path)

    # make directory for plots
    print('making dir for figures')
    figures_path = results_dir + 'Figures/'
    try:
        makedirs(figures_path)
    except FileExistsError:
        pass

    # generate plots
    print('generating reference and region plot')
    ref_bin_counts = plot_59_genes_reference_gc_bias(bed, ref, figures_path)
    print('generating plots for samples')
    make_plots(summary_stats_path, figures_path, ref_bin_counts)


"""
Running from command line parameters
1. reference genome
2. path directory containing all bam files
3. path directory containing all bed
4. path of the results directory
"""


if __name__ == "__main__":
    # get all the bams in the bam dir
    ref = sys.argv[1]
    bam_template = sys.argv[2]
    bed_template = sys.argv[3]
    res = sys.argv[4]

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
    run_gc_bias_analysis(ref, bams, beds, res)
