import argparse
import os
import numpy as np

#sample_qualities_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/metric_files/'
#hm_metrics_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/heatmap_metrics/'

def write_quality_per_gene():
    # store arguments from command line
    args = get_cmdln_arguments()
    sample_metrics_path = args['m'][0] 
    hm_metrics_path = args['hm'][0]

    for sample in os.listdir(sample_metrics_path):
        print(sample)
        # read in all qualities
        curr_sample = open(sample_metrics_path + sample, 'r')

        # get avg quality per gene
        gene_qualities = get_avg_quality_per_gene(curr_sample)

        # write avg quality for each gene
        heatmap_output = open(hm_metrics_path + sample.split('.')[0] + '_hm.txt', 'a')

        for gene in gene_qualities.keys():
            heatmap_output.write(gene + '\t' + str(gene_qualities[gene]) + '\n')

def get_avg_quality_per_gene(curr_sample):
    sample_qualities_dict = {}
    for line in curr_sample:
        line = line.rstrip().split()
        try: sample_qualities_dict[line[4]].append(float(line[3]))
        except KeyError: sample_qualities_dict[line[4]] = [float(line[3])]
    for gene in sample_qualities_dict.keys():
        sample_qualities_dict[gene] = np.average(sample_qualities_dict[gene])
    print(sample_qualities_dict)
    return sample_qualities_dict

# get the arguments from commandline run
def get_cmdln_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', nargs = 1, required = True, help = 'specify the directory of metric files')
    parser.add_argument('-hm', nargs = 1, required = True, help = 'specify the directory of heatmap files (output)')
    args = parser.parse_args()
    return vars(args)


write_quality_per_gene()
