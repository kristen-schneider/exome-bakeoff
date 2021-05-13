import argparse
import glob
import os
import collections
import matplotlib.pyplot as plt
import numpy as np


# input text file with all directories of final metrics for each sample
# sample name
def get_args():

    parser = argparse.ArgumentParser()


    parser.add_argument('--in_txt',
                        dest='in_txt',
                        type=str,
                        help='path to txt file with all directories storing mean_per_tech_combo results',
                        required=True)
    parser.add_argument('--sample',
                        dest='sample',
                        help='sample name wanted to compute over (e.g. GIAB, etc.). Default is all.',
                        required=False)
        
    parser.add_argument('--out_csv',
                        dest='out_csv',
                        help='Path of output csv file',
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
    # list of all metrics
    metric_names = get_metric_names(args.in_txt)
    #tech_combo_dict = get_metric_dict(args.in_txt, args.sample)
    #print(tech_combo_dict)
    raw_scores_dict = raw_metric_dict(args.in_txt, args.sample)
    print(raw_scores_dict)
    write_out_file(raw_scores_dict, args.in_txt, args.out_csv)
    #plot_final_heatmap(tech_combo_dict, args.in_txt, args.out_png)

# takes a list of directories for all metrics and returns a list of the metrics evaluated
def get_metric_names(in_txt):
    metrics_names = []
    for line in open(in_txt):
        metric = line.rstrip().split('/')[8]
        if metric == 'F1score': metric = "F1_"+line.rstrip().split('/')[9].split('_')[1]
        metrics_names.append(metric)
    
    return metrics_names

def transpose_dict_data(tech_combo_dict):
    data = []
    for i in tech_combo_dict:
        data.append(tech_combo_dict[i])
    
    numpy_array = np.array(data)
    transpose = numpy_array.T
    transpose_list = transpose.tolist()
    return transpose_list

def plot_final_heatmap(tech_combo_dict, in_txt, out_png):
    tech_combos_names = []
    for i in tech_combo_dict: tech_combos_names.append(i)
    metrics_names = get_metric_names(in_txt)

    data = transpose_dict_data(tech_combo_dict)

    plt.figure(figsize=(20, 10))
    plt.imshow(data, interpolation='nearest', cmap='Greens', vmin=0, vmax=1)#, interpolation=None, cmap='Blues')#, vmin=0, vmax=1)

    plt.colorbar()
    plt.title('Metrics for 59 ACMG Genes for GIAB samples')
    plt.yticks(np.arange(0, len(metrics_names)), metrics_names, rotation=0)
    plt.ylabel('METRICS')
    plt.xticks(np.arange(0, len(tech_combos_names)), tech_combos_names, rotation=90)
    plt.xlabel('TECH-SAMPLE COMBO')
    plt.savefig(out_png)


def write_out_file(tech_combo_dict, in_txt, out_csv):
    metrics_names = get_metric_names(in_txt)
    
    f = open(out_csv, 'w')
    f.write('')
    for m in metrics_names: f.write(','+m)
    f.write('\n')
    for i in tech_combo_dict:
        f.write(i+','+str(tech_combo_dict[i]).replace('[', '').replace(']','')+'\n')

def get_max_score(metric_dir, sample, curr_metric):
    if curr_metric == 'quality':
        print('quality')
        max_score = simple_max(metric_dir, sample)
    elif curr_metric == 'strandbias':
        print('strandbias')
        max_score = difference_max(metric_dir, sample, 0.0)
    elif curr_metric == 'noise':
        print('noise')
        max_score = simple_max(metric_dir, sample)
    elif curr_metric == 'gcbias':
        print('gcbias')
        max_score = difference_max(metric_dir, sample, 0.0)
    elif curr_metric == 'F1score':
        print('F1score')
        max_score = simple_max(metric_dir, sample)

    else: max_score = None 
    return max_score

# returns max of list
# where list is all scores for all samples of a given metric_dir
def simple_max(metric_dir, sample):
    all_metric_scores = []
    for tech_combo_file in os.listdir(metric_dir):
        if sample in tech_combo_file and ".txt" in tech_combo_file:
            for line in open(metric_dir+tech_combo_file):
                score = float(line.rstrip())
                all_metric_scores.append(score)
    return max(all_metric_scores)

# finds difference between score and some base value,
# then returns max of those scores
def difference_max(metric_dir, sample, base_value):
    all_metric_scores = []
    for tech_combo_file in os.listdir(metric_dir):
        if sample in tech_combo_file and ".txt" in tech_combo_file:
            for line in open(metric_dir+tech_combo_file):
                score = abs(base_value - float(line.rstrip()))
                all_metric_scores.append(score)
    return min(all_metric_scores)

def raw_metric_dict(in_txt, sample):
    tech_combo_dict = dict()

    # process specified samples
    if sample == None: token = '.txt'
    else: token = sample

    # iterate through one metric directory
    for line in open(in_txt):
        metric_dir = line.rstrip()
        curr_metric = metric_dir.split('/')[8]
        print(metric_dir)
    
        for tech_combo_file in glob.glob(metric_dir+'/*'):
            if token in tech_combo_file and '.txt' in tech_combo_file:
                basename = os.path.basename(tech_combo_file)
                tech_combo_name = basename.split('_')[0]
                if tech_combo_name not in tech_combo_dict: tech_combo_dict[tech_combo_name] = []

                for line in open(tech_combo_file):
                    score = abs(float(line.rstrip()))
                    #print(tech_combo_name+': '+str(score))
                    tech_combo_dict[tech_combo_name].append(score)
    sorted_dict = collections.OrderedDict(sorted(tech_combo_dict.items()))
    return sorted_dict    


# returns a sorted dictinary {key: sample name, value: [list of scores for a metric]}
def get_metric_dict(in_txt, sample):
    tech_combo_dict = dict()
    
    # process specified samples
    if sample == None: token = '.txt'
    else: token = sample 

    # iterate through one metric directory
    for line in open(in_txt):
        metric_dir = line.rstrip()
        curr_metric = metric_dir.split('/')[8]
        print(metric_dir)
        
        # find max value in this sample set (to use for normalize)
        max_score = get_max_score(metric_dir, sample, curr_metric)
        print(max_score)    
        for tech_combo_file in glob.glob(metric_dir+'/*'):
            if token in tech_combo_file and '.txt' in tech_combo_file:
                basename = os.path.basename(tech_combo_file)
                tech_combo_name = basename.split('_')[0]
                if tech_combo_name not in tech_combo_dict: tech_combo_dict[tech_combo_name] = []
            
                for line in open(tech_combo_file):
                    score = abs(float(line.rstrip()))
                    print(tech_combo_name+': '+str(score))
                    # reported_score / max_score for this metric
                    if curr_metric == 'quality' or curr_metric == 'noise' or curr_metric == 'F1score':
                        tech_combo_dict[tech_combo_name].append(score/max_score)
                    else:
                        tech_combo_dict[tech_combo_name].append(max_score/score)
            else: continue
    
    sorted_dict = collections.OrderedDict(sorted(tech_combo_dict.items()))
    return sorted_dict
                



if __name__ == '__main__':
    main()

