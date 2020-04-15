import os
import numpy as np
import matplotlib.pyplot as plt

sample_qualities_path = '/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/metric_files/'
scatter_plot_path = '/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/scatterplots/'
TITLE = 'Quality'

def main_scatter():
    for sample in os.listdir(sample_qualities_path):
        # read in all qualities
        curr_sample = open(sample_qualities_path + sample, 'r')

        # get qualities per gene
        gene_qualities = get_qualities_per_gene(curr_sample)

        # plot one scatter for each gene
        for gene in gene_qualities.keys():
            plot_scatter(gene, gene_qualities[gene], sample)

def get_qualities_per_gene(curr_sample):
    sample_qualities_dict = {}
    for line in curr_sample:
        line = line.rstrip().split()
        pos_qual = [int(line[1]), float(line[3])]
        try: sample_qualities_dict[line[4]].append(pos_qual)
        except KeyError: sample_qualities_dict[line[4]] = [pos_qual]
    return sample_qualities_dict

def plot_scatter(gene, gene_qualities, sample):
    sample = sample.split('.')[0]
    x = []
    for g in gene_qualities: x.append(g[0])
    y = []
    for g in gene_qualities: y.append(g[1])

    plt.figure(figsize=(35, 10))
    plt.plot(x, y, 'o', color='black')
    plt.ylim((0,93))
    y_label = 'Quality'


    plt.title(y_label + ' for gene ' + gene + ' in sample ' + sample)
    plt.xlabel('Exome Position')
    plt.ylabel(y_label)
    plt.savefig(scatter_plot_path + sample+'-'+gene+'scatter.png')

main_scatter()