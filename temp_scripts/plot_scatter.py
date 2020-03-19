import sys
import os
import numpy as np
import matplotlib.pyplot as plt

final_metrics = sys.argv[1]
gene = sys.argv[2]
sample = sys.argv[3]

def main():
    all_data = []
    for metric_file in os.listdir(final_metrics):
        if gene in metric_file and sample in metric_file:
            data_points = read_file(final_metrics, metric_file)
            scatter_plot(metric_file, data_points, gene, sample)

def read_file(final_metric, metric_file):
    f = open(os.path.join(final_metrics, metric_file), 'r')
    x = []
    y = []
    #full_xy = []
    for line in f:
        line = line.strip().split()
        x.append(float(line[0].strip()))
        y.append(float(line[1].strip()))
    full_xy = [x, y]
    return full_xy

def scatter_plot(metric_file, data_points, gene, sample):
    x = data_points[0]
    y = data_points[1]
    
    plt.figure(figsize=(35, 10))
    plt.plot(x, y, 'o', color='black')
    y_label = ''
    if 'quality' in metric_file:
        plt.ylim((0,93))
        y_label = 'Quality'
    elif '_bias' in metric_file:
        y_label = 'Strand Bias'
        plt.axhline(y=get_forward_strand_mean_std(data_points), color='red')
        plt.axhline(y=get_reverse_strand_mean_std(data_points), color='red')
    
    plt.title(y_label + ' at position x for gene ' + gene + ' in sample ' + sample) 
    plt.xlabel('Position')
    plt.ylabel(y_label)
    plt.savefig(os.path.basename(y_label+sample+'-'+gene+'scatter.png'))

def get_forward_strand_mean_std(xy_list):
    y = xy_list[1]
    forwards = []
    for i in y:
        if i < 0.5: forwards.append(i)
    return np.average(forwards)-2*np.std(forwards)

def get_reverse_strand_mean_std(xy_list):
    y = xy_list[1]
    reverses = []
    for i in y:
        if i >= 0.5: reverses.append(i)
    return np.average(reverses)+2*np.std(reverses)

main()
