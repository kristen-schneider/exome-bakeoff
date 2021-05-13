import argparse
import os
import operator
import matplotlib.pyplot as plt
import numpy as np


def get_args():

    parser = argparse.ArgumentParser()


    parser.add_argument('--in_csv',
                        dest='in_csv',
                        type=str,
                        help='path to csv file with raw scores for each metric',    
                        required=True)
    
    parser.add_argument('--out_csv',
                        dest='out_csv',
                        help='Path of output csv ranking file',
                        type=str,
                        required=True)

    parser.add_argument('--out_png',
                        dest='out_png',
                        help='Path of output png heatmap file',
                        type=str,
                        required=True)

    return parser.parse_args()

def main():
    args = get_args()
    metrics = get_metric_names(args.in_csv)
    csv_rank = open(args.out_csv, 'w')
    for m in metrics: csv_rank.write(str(m)+',')
    csv_rank.write('\n')
    
    for m in metrics: 
        if m != '':
            if m == 'quality':
                quality = make_dict(args.in_csv, metrics.index(m)) 
                quality_rank = rank_dict(quality)
            elif m == 'strandbias':
                sb = make_dict(args.in_csv, metrics.index(m))
                sb_rank = rank_dict(sb)
            elif m == 'noise':
                noise = make_dict(args.in_csv, metrics.index(m))
                noise_rank = rank_dict(noise)
            elif m == 'gcbias':
                gcbias = make_dict(args.in_csv, metrics.index(m))
                gc_bias_rank = rank_dict(gcbias)
            elif m == 'F1_INDEL':
                indel = make_dict(args.in_csv, metrics.index(m))
                indel_rank = rank_dict(indel)
            elif m == 'F1_SNP':
                snp = make_dict(args.in_csv, metrics.index(m))
                snp_rank = rank_dict(snp)                
    for t in indel:
        csv_rank.write(t + ',' + str(indel_rank[t]) + ',' + str(snp_rank[t])+'\n')
    
    #for t in quality:
    #    csv_rank.write(t + ',' + str(quality_rank[t]) + ',' + str(sb_rank[t]) + ',' + str(noise_rank[t]) + ',' + str(gc_bias_rank[t])+'\n')
    #print(metric_dict)
    #print(metric_dict_rank)            

def rank_dict(metric_dict):
    sorted_dict = dict()
    values = metric_dict.values()
    sorted_values = sorted(values, reverse=True)
    for t in metric_dict:
        for s in sorted_values:
            if metric_dict[t] == s:
                sorted_dict[t] = sorted_values.index(s)+1
    #print(values)
    #print(sorted_values)
    return sorted_dict


def get_metric_names(in_csv):
    header = ''
    for line in open(in_csv):
        if header == '':
            header = line.rstrip().split(',')
        break
    return header

# makes a dictionary for one metric (sample: score)
def make_dict(in_csv, index):
    header = ''
    metric_dict = dict()
    for line in open(in_csv):
        if header == '': header = line
        else:
            A = line.rstrip().split(',')
            tech_combo = A[0]
            metric_score = float(A[index])
            metric_dict[tech_combo] = metric_score
    return metric_dict

# https://www.geeksforgeeks.org/python-program-to-sort-a-list-of-tuples-by-second-item/
def Sort_Tuple(tup):

    # getting length of list of tuples 
    lst = len(tup)
    for i in range(0, lst):

        for j in range(0, lst-i-1):
            if (tup[j][1] > tup[j + 1][1]):
                temp = tup[j]
                tup[j]= tup[j + 1]
                tup[j + 1]= temp
    return tup


def plot_ranks(metrics, tech_combos, data, out_png):
    metrics = metrics[1:]
    plt.figure(figsize=(20, 10))
    plt.imshow(data, interpolation='nearest', cmap='Greens', vmin=0, vmax=len(tech_combos)-1)#, interpolation=None, cmap='Blues')#, vmin=0, vmax=1)

    plt.colorbar()
    plt.title('Metrics for 59 ACMG Genes for GIAB samples')
    plt.yticks(np.arange(0, len(metrics)), metrics, rotation=0)
    plt.ylabel('METRICS')
    plt.xticks(np.arange(0, len(tech_combos)), tech_combos, rotation=90)
    plt.xlabel('TECH-SAMPLE COMBO')
    plt.savefig(out_png)

    return(ranks_by_tech_combo) 
 
    
def get_tech_combos(in_csv):
    tech_combos = []
    csv = open(in_csv, 'r')
    header = ''
    for line in csv:
        if header == '': header = line
        else: tech_combos.append(line.rstrip().split(',')[0])
    return tech_combos
         

if __name__ == '__main__':
    main()
