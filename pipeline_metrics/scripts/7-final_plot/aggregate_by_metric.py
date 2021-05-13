import argparse
import glob
import os
import collections
import matplotlib.pyplot as plt
import numpy as np

def get_args():

    parser = argparse.ArgumentParser()


    parser.add_argument('--in_txt',
                        dest='in_txt',
                        type=str,
                        help='path to txt file with all directories storing mean_per_tech_combo results',
                        required=True)

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
    tech_combo_dict = get_dict(args.in_txt)
    write_out_file(tech_combo_dict, args.out_csv)
    plot_final_heatmap(tech_combo_dict, args.in_txt, args.out_png)

def get_metrics_names(in_txt):
    metrics_names = []
    for line in open(in_txt):
        metric = line.rstrip().split('/')[8]
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
    metrics_names = get_metrics_names(in_txt)

    data = transpose_dict_data(tech_combo_dict)

    plt.figure(figsize=(40, 5))
    plt.imshow(data, interpolation='nearest', cmap='Greens', vmin=0, vmax=1)#, interpolation=None, cmap='Blues')#, vmin=0, vmax=1)

    plt.colorbar()
    plt.yticks(np.arange(0, len(metrics_names)), metrics_names, rotation=0)
    plt.ylabel('METRICS')
    plt.xticks(np.arange(0, len(tech_combos_names)), tech_combos_names, rotation=90)
    plt.xlabel('TECH-SAMPLE COMBO')
    plt.savefig(out_png)


def write_out_file(tech_combo_dict, out_csv):
    f = open(out_csv, 'w')
    for i in tech_combo_dict:
        f.write(i+','+str(tech_combo_dict[i]).replace('[', '').replace(']','')+'\n')

def get_max_score(metric_dir):
    all_metric_scores = []
    for tech_combo_file in os.listdir(metric_dir):
        if ".txt" in tech_combo_file:
            for line in open(metric_dir+tech_combo_file):
                score = float(line.rstrip())
                all_metric_scores.append(score)
    #print(all_metric_scores)
    return max(all_metric_scores)
        
 
def get_dict(in_txt):
    tech_combo_dict = dict()
    
    for line in open(in_txt):
        metric_dir = line.rstrip()
        print(metric_dir)
        max_score = get_max_score(metric_dir)

        for tech_combo_file in glob.glob(metric_dir+'/*'):
            if ".txt" in tech_combo_file:
                basename = os.path.basename(tech_combo_file)
                tech_combo_name = basename.split('_')[0]
                if tech_combo_name not in tech_combo_dict: tech_combo_dict[tech_combo_name] = []
            
                for line in open(tech_combo_file):
                    score = float(line.rstrip())
                    # reported_score / max_score for this metric
                    tech_combo_dict[tech_combo_name].append(score/max_score)
            else: continue
    
    sorted_dict = collections.OrderedDict(sorted(tech_combo_dict.items()))
    

    return sorted_dict
                



if __name__ == '__main__':
    main()

