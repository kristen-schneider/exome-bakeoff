import argparse
import os
import numpy as np
import matplotlib.pyplot as plt


#heatmap_metrics_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/heatmap_metrics/'
#heatmaps_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/'
#TITLE = 'Quality'

def main_heatmap():
    # store arguments from command line
    args = get_cmdln_arguments()
    heatmap_metrics_path = args['hm'][0]
    heatmaps_figure_path = args['o'][0]
    TITLE = args['t'][0]

    sample_names = []
    gene_names = []
    all_sample_metrcis = []

    for sample in os.listdir(heatmap_metrics_path):
        if "_hm" in sample:
            # store sample names
            curr_sample_name = sample.split('-')
            curr_sample_name = curr_sample_name[0]+'-'+curr_sample_name[1]+'-'+curr_sample_name[2]
            sample_names.append(curr_sample_name)

            # store gene names with data
            curr_sample_gene_data = {}
            curr_sample_file = open(heatmap_metrics_path + sample, 'r')
            for line in curr_sample_file:
                line = line.rstrip().split()
                gene = line[0]
                metric = line [1]
                if gene not in gene_names: gene_names.append(gene)
                curr_sample_gene_data[gene] = metric
            # sort sample metrics by gene name
            genes_sorted = sorted(curr_sample_gene_data)
            #print(genes_sorted)

            metrics_sorted_by_gene = []
            for gene in genes_sorted: metrics_sorted_by_gene.append(float(curr_sample_gene_data[gene]))

            all_sample_metrcis.append(metrics_sorted_by_gene)

    max = 0
    min = 100
    for m in all_sample_metrcis:
        for i in m:
            if i > max: max = i
            if i < min: min = i
    print(max)
    print(min)

    plot_heatmap(sample_names, gene_names, all_sample_metrcis, TITLE, heatmaps_figure_path)

def plot_heatmap(sample_names, gene_names, all_sample_metrcis, title, heatmaps_figure_path):
    fig, ax = plt.subplots()
    fig.figsize = (40,30)
    fig.set_size_inches(35.0, 25.0)
    ax.set_xticks(np.arange(len(gene_names)))
    ax.set_yticks(np.arange(len(sample_names)))
    ax.set_xticklabels(gene_names)
    ax.set_yticklabels(sample_names)
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.imshow(all_sample_metrcis, interpolation=None, cmap='coolwarm', vmin=35, vmax=58)
    plt.colorbar()
    plt.savefig(heatmaps_figure_path + title+'.png', dpi=150)

# get the arguments from commandline run
def get_cmdln_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hm', nargs = 1, required = True, help = 'specify the directory of heatmap files')
    parser.add_argument('-o', nargs = 1, required = True, help = 'specify the directory of figure output')
    parser.add_argument('-t', nargs = 1, required = True, help = 'give the plot a title')
    args = parser.parse_args()
    return vars(args)

main_heatmap()
