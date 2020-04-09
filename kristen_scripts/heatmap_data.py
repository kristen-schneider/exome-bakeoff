import os
import numpy as np

sample_qualities_path = '/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/metric_files/'
hm_metrics_path = '/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/heatmap_metrics/'

def write_quality_per_gene():
    for sample in os.listdir(sample_qualities_path):
        # read in all qualities
        curr_sample = open(sample_qualities_path + sample, 'r')

        # get qualities per gene
        gene_qualities = get_quality_per_gene(curr_sample)

        # write avg quality for each gene
        heatmap_output = open(hm_metrics_path + sample.split('.')[0] + '-hm.txt', 'a')

        for gene in gene_qualities.keys():
            heatmap_output.write(gene + '\t' + str(gene_qualities[gene]) + '\n')

def get_quality_per_gene(curr_sample):
    sample_qualities_dict = {}
    for line in curr_sample:
        line = line.rstrip().split()
        try: sample_qualities_dict[line[4]].append(float(line[3]))
        except KeyError: sample_qualities_dict[line[4]] = [float(line[3])]
    for gene in sample_qualities_dict.keys():
        sample_qualities_dict[gene] = np.average(sample_qualities_dict[gene])
    return sample_qualities_dict


write_quality_per_gene()
