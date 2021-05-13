import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import collections
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


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
                        type=str)#,
                        #required=True)

    parser.add_argument('--name',
                        dest='name',
                        type=str)#,
                        #required=True)    

    return parser.parse_args()


def main():
    args = get_args()
    preps = get_preps(args.input_path)
    captures = get_captures(args.input_path)
    hm_dict = make_hm_dict(args.input_path)   
    hm_data = make_hm_data(preps, captures, hm_dict)
    plot_hm(preps, captures, hm_data)
    
    
def plot_hm(preps, captures, hm_data): 
    plt.figure(figsize=(8, 6))
    
    # ryan's code to make null elements greyed out
    viridis = cm.get_cmap('Greens')
    newcolors = viridis(np.linspace(0, 1, 256))
    grey = np.array([0.5, 0.5, 0.5, 0.5])
    newcolors[:25, :] = grey
    newcmp = ListedColormap(newcolors)
    #
    
    plt.imshow(hm_data, interpolation='nearest', cmap=newcmp ,vmin=0.3, vmax=0.5)#, interpolation=None, cmap='Blues')#, vmin=0, vmax=1)
        
    
    # annotating heatmap
    
    for i in range(len(hm_data)):
        for j in range(len(hm_data[i])):
            #print(data[i][j])
            if hm_data[i][j] > 0: text = plt.text(j, i, round(hm_data[i][j], 3),
                       ha="center", va="center", color="w")

    
    plt.colorbar()
    plt.xticks(np.arange(0, len(preps)), preps, rotation=65)
    plt.xlabel('LIBRARY PREP')
    plt.yticks(np.arange(0, len(captures)), captures, rotation=65)
    plt.ylabel('CAPTURE')
    plt.savefig('tech_hm.png')


def make_hm_data(preps, captures, hm_dict):
    empty = -1
    hm_data = [[empty]*len(preps) for n in range(len(captures))]
    
    for tech_combo in hm_dict:
        t_prep = tech_combo.split('-')[0]
        t_capture = tech_combo.split('-')[1]
        
        hm_data[captures.index(t_capture)][preps.index(t_prep)] = hm_dict[tech_combo]
    return hm_data   

def make_hm_dict(input_path):
    tech_sample_dict = dict()
    
    for f in os.listdir(input_path): 
        tech_sample = ('-').join(f.split('_')[0].split('-')[0:2])
        if tech_sample not in tech_sample_dict: tech_sample_dict[tech_sample] = [] 
        o = open(input_path + f, 'r')
        
        tech_sample_mean = 0
        tech_sample_scores = []
        
        for line in o:
            A = line.rstrip().split()
            gene = A[0]
            score = float(A[1])
            tech_sample_scores.append(score)
        tech_sample_mean = sum(tech_sample_scores)/len(tech_sample_scores)
        tech_sample_dict[tech_sample] = tech_sample_mean
    return tech_sample_dict  

 
def get_preps(input_path):
    preps = []
    for f in os.listdir(input_path):
        prep = f.split('-')[0]
        if prep not in preps: preps.append(prep)
    return sorted(preps)

def get_captures(input_path):
    captures = []
    for f in os.listdir(input_path):
        capture = f.split('-')[1]
        if capture not in captures: captures.append(capture)
    return sorted(captures)





if __name__ == '__main__':
    main()
