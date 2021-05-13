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
                        help='path to csv file with rankings for each metric',
                        required=True)

    parser.add_argument('--out_png',
                        dest='out_png',
                        help='Path of output png heatmap file',
                        type=str,
                        required=True)

    return parser.parse_args()

def main():
    args = get_args()
    data = read_rank_csv(args.in_csv)
    plot_ranks(data, args.out_png)

def read_rank_csv(rank_csv):
    metrics = []
    tech_combos = []
    data = []
    for line in open(rank_csv):
        if len(metrics) == 0: metrics = line.rstrip().split(',')[1:]
        else:
            A = line.rstrip().split(',')
            tech = A[0]
            r = A[1:]
            rankings = []
            for i in r: rankings.append(int(i))
            tech_combos.append(tech)
            data.append(rankings)
    transposed = transpose_data(data)
    # print(metrics)
    # print(tech_combos)
    # print(transposed)
    return [metrics, tech_combos, transposed]

def transpose_data(data):
    transposed = []
    for m in data[0]:
        transposed.append([])

    for t in data:
        for m in range(len(t)):
            metric_index = m
            value_to_append = t[m]
            transposed[metric_index].append(value_to_append)
    return transposed

def plot_ranks(data, out_png):
    metrics = data[0][0:2]
    print(metrics)
    tech_combos = data[1]
    plot_data = data[2]

    plt.figure(figsize=(20, 10))
    plt.imshow(plot_data, interpolation='nearest', cmap='Greens', vmin=1, vmax=len(tech_combos))#, interpolation=None, cmap='Blues')#, vmin=0, vmax=1)

    plt.colorbar()
    plt.title('F1 Score for GIAB samples')
    plt.yticks(np.arange(0, len(metrics)), metrics, rotation=0)
    plt.ylabel('METRICS')
    plt.xticks(np.arange(0, len(tech_combos)), tech_combos, rotation=90)
    plt.xlabel('TECH-SAMPLE COMBO')
    plt.savefig(out_png)

    #return(ranks_by_tech_combo)

if __name__ == '__main__':
    main()
