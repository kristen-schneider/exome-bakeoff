import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re
import sys
from os import makedirs


def plot_by_tech(df, col, y_col, y_lab, fig_name):
    # rgba based color schemes to make colors semi transparent
    # pale = [(.07,.46,.2,.5), (.2,.14,.53,.5),(.86,.8,.46,.5),(.8,.4,.46,.5),(.53,.01,.33,.5),(.54,.8,.93,.5)]
    # opacity = .8
    # colors = [(230/255,159/255,3/255,opacity),(86/255,180/255,233/255,opacity),(5/255,158/255,115/255,opacity),(3/255,114/255,178/255,opacity),(213/255,94/255,0,opacity),(204/255,121/255,167/255,opacity)]
    colors = ['r', 'g', 'b', 'orange', 'purple', 'y']
    widths = [6, 5, 4, 3, 2, 1]
    samples = list(set(df[col]))
    # create dictionaries of the group and the color / width of their lines
    color_map = {samples[i]: colors[i] for i in range(len(samples))}
    width_map = {samples[i]: widths[i] for i in range(len(samples))}
    custom_lines = []

    ax = plt.subplot(111)
    for sample in samples:
        custom_lines.append(Line2D([0], [0], color=color_map[sample], lw=4))
        sub = df[df[col] == sample]
        for s in list(set(sub['SAMPLE'])):
            subsub = sub[sub['SAMPLE'] == s]
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


def make_plots(csv_dir, results_dir):
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

    # create alternate columns to plot with
    data['ref_bin_normed_count'] = data['count'] / data['ref_windows']
    data['tech_one'] = [re.sub('-.*-.*', '', x) for x in list(data['SAMPLE'])]
    data['tech_two'] = [re.sub('.*-(\\w+)-.*', '\\1', x) for x in list(data['SAMPLE'])]

    # line plot with every sample as its own color
    ax = sns.lineplot(x="percent", y="ref_bin_normed_count", data=data, hue='SAMPLE')
    ax.legend_.remove()
    plt.savefig(results_dir + 'ref_bin_normalized.png')

    # ribbon plot with ribbon representing standard deviation grouped by library prep
    ax = sns.lineplot(x="percent", y="ref_bin_normed_count", data=data, hue='tech_one', ci='sd')
    ax.legend_.remove()
    plt.savefig(results_dir + 'ribbon-gc_bias-library_prep.png')

    # ribbon plot with ribbon representing standard deviation grouped by capture technology
    ax = sns.lineplot(x="percent", y="ref_bin_normed_count", data=data, hue='tech_two', ci='sd')
    ax.legend_.remove()
    plt.savefig(results_dir + 'ribbon-gc_bias-capture_technology.png')
    plt.clf()

    plot_by_tech(data, 'tech_one', 'ref_bin_normed_count', 'Normalized Read Count', results_dir + 'gc_bias_library_prep.png')
    plot_by_tech(data, 'tech_two', 'ref_bin_normed_count', 'Normalized Read Count', results_dir + 'gc_bias_capture_technology.png')

    plot_by_tech(data, 'tech_one', 'count', 'Read Count', results_dir + 'raw_coverage_gc_bias_library_prep.png')
    plot_by_tech(data, 'tech_two', 'count', 'Read Count', results_dir + 'raw_coverage_gc_bias_capture_technology.png')


if __name__ == "__main__":
    input_dir = sys.argv[1]
    res_dir = sys.argv[2]
    make_plots(input_dir, res_dir)
