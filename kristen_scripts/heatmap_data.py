import os
import numpy as np

sample_qualities_path = '/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/metric_files/'

def write_quality_per_gene():
    for sample in os.listdir(sample_qualities_path):
        curr_sample = open(sample_qualities_path + sample, 'r')

        heatmap_output = open(sample+curr_sample)

def get_quality_per_gene(curr_sample):
    sample_qualities_dict = {}
    for line in curr_sample:
        line = line.rstrip().split()
        try: sample_qualities_dict[line[4]].append(float(line[3]))
        except KeyError: sample_qualities_dict[line[4]] = [float(line[3])]
    for gene in sample_qualities_dict.keys():
        print(gene)
        print(np.average(sample_qualities_dict[gene]))
        sample_qualities_dict[gene] = np.average(sample_qualities_dict[gene])
    return sample_qualities_dict


