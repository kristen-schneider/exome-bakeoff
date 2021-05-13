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

    parser.add_argument('--tech',
                        dest='tech',
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
    if (args.tech == 'prep'):
        preps_data = get_preps_data(args.input_path)
        plot_boxplot(preps_data)
        check(preps_data)
    if (args.tech == 'capture'):
        capture_data = get_capture_data(args.input_path)
        plot_boxplot(capture_data)
        check(capture_data)

def get_preps_data(input_path):
    # holds all values for all genes in the given library prep
    prep_dict = dict()
    for tech_sample_combo in os.listdir(input_path):
        prep_name = tech_sample_combo.split("-")[0]
        f = open(input_path+tech_sample_combo, 'r')
        for line in f:
            A = line.rstrip().split()
            gene = A[0]
            score = float(A[1])
            try: prep_dict[prep_name].append(score)
            except KeyError: prep_dict[prep_name] = [score]
    return prep_dict            

def get_capture_data(input_path):
    # holds all values for all genes in the given capture
    capture_dict = dict()
    for tech_sample_combo in os.listdir(input_path):
        capture_name = tech_sample_combo.split("-")[1]
        f = open(input_path+tech_sample_combo, 'r')
        for line in f:
            A = line.rstrip().split()
            gene = A[0]
            score = float(A[1])
            try: capture_dict[capture_name].append(score)
            except KeyError: capture_dict[capture_name] = [score]
    return capture_dict

def check(data):
    for t in data.keys():
        print(t)
        print(sum(data[t])/len(data[t]))
        print(max(data[t]))
        print(min(data[t]))


def plot_boxplot(data):
    # sort alphabetically and get techs      
    od = collections.OrderedDict(sorted(data.items()))
    tech_list = []
    data_by_tech = []
    for t in od: tech_list.append(t)
   
    techs = []
    data = []
    for t in od:
        techs.append(t)
        data.append(od[t])
    data_by_tech = [techs, data]

    plt.figure(figsize=(15, 10))
    box = plt.boxplot(data_by_tech[1])
    plt.xticks(np.arange(1, len(data_by_tech[0])+1), data_by_tech[0], rotation=65)
    #plt.xticks(np.arange(1, len(tech_list)), tech_list, rotation=65) 

    plt.savefig('boxplot.png')
    

if __name__ == '__main__':
    main()
