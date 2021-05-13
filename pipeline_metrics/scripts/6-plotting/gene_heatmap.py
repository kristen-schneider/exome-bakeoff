import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import collections

def get_args():
   
    # input the path to the metric files
    # format: gene  metric
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path',
                        dest='input_path',
                        type=str,
                        required=True)

    parser.add_argument('--output_path',
                        dest='output_path',
                        type=str,
                        required=True)

    parser.add_argument('--name',
                        dest='name',
                        type=str,
                        required=True)    

    return parser.parse_args()


def main():
    args = get_args()
    sorted_tech_genes = get_sorted_tech_genes(args.input_path)
    tech_sample_gene_list = make_tech_sample_gene_list(args.input_path)
    gene_heatmap_data = make_gene_heatmap_data(tech_sample_gene_list)
    plot_gene_heatmap(gene_heatmap_data, sorted_tech_genes, args.output_path, args.name)#, basename, out_dir) 
    
def make_tech_sample_gene_list(input_path):
    tech_sample_list = []

    # go through all files (tech_sample combos) in dir
    for metric_file in os.listdir(input_path):
        #print(metric_file)  
        tech_sample_combo = metric_file.split('_')[0]
        
        # go through all lines (genes) in file
        f = open(input_path+metric_file, 'r')
        gene_score_dict = dict()
        for line in f:
            A = line.rstrip().split()
            gene = A[0]
            score = float(A[1])
            gene_score_dict[gene] = score
            #gene_score_list.append((gene, score))
            od = collections.OrderedDict(sorted(gene_score_dict.items()))       

     
        tech_sample_list.append((tech_sample_combo, od))
    return tech_sample_list
    
def make_gene_heatmap_data(tech_sample_list):
    data = []
    
    sorted_by_tech = sorted(tech_sample_list)
    for tech_sample_combo in sorted_by_tech:
        gene_scores = []
        for gene_score in tech_sample_combo[1]:
            gene_scores.append(tech_sample_combo[1][gene_score])
            #print(tech_sample_combo[1][gene_score])
        #for gene in tech_sample_combo:
            #print(gene)
            #gene_scores.append(tech_sample_combo[1][gene])
        
        data.append(gene_scores)
    return data

# make a list of tech-samples and genes in sorted order
def get_sorted_tech_genes(input_path):
    tech_sample_list = []
    gene_list = []
    # go through all files (tech_sample combos) in dir
    for metric_file in os.listdir(input_path):
        #print(metric_file)  
        tech_sample_combo = metric_file.split('_')[0]
        if tech_sample_combo not in tech_sample_list: tech_sample_list.append(tech_sample_combo)
        # go through all lines (genes) in file
        f = open(input_path+metric_file, 'r')
        for line in f:
            A = line.rstrip().split()
            gene = A[0]
            if gene not in gene_list: gene_list.append(gene)
    return [sorted(tech_sample_list), sorted(gene_list)]

def plot_gene_heatmap(data, sorted_tech_genes, out_dir, name):#, basename, out_dir):
    plt.figure(figsize=(40, 50))
    plt.imshow(data, interpolation='nearest', cmap='Greens', vmin=0, vmax=1)#, interpolation=None, cmap='Blues')#, vmin=0, vmax=1)
        
    
    plt.colorbar()
    plt.xticks(np.arange(0, len(sorted_tech_genes[1])), sorted_tech_genes[1], rotation=65)
    plt.xlabel('GENES')
    plt.yticks(np.arange(0, len(sorted_tech_genes[0])), sorted_tech_genes[0], rotation=0)
    plt.ylabel('TECH-SAMPLE COMBO')
    plt.savefig(out_dir + name+'.png')



    


if __name__ == '__main__':
    main()
