import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def pm_main(final_metrics_path, region_file, num_samples):
    curr_gene = os.path.splitext(region_file)[0]
    
    
    for txt_file in os.listdir(final_metrics_path):
        if curr_gene in txt_file and txt_file.endswith('final.txt'):
            xy_list = read_file(final_metrics_path, txt_file)
            #scatter_plot(xy_list, txt_file)
            if txt_file.endswith('quality_final.txt'): write_heat_map_avg_residual(txt_file, region_file, xy_list)
            if txt_file.endswith('strand_bias_final.txt'): write_heat_map_proportion_outside(txt_file, region_file, xy_list) 
    #plot_heat_map() 
    
def read_file(final_metrics_path, metrics_txt):
    f = open(os.path.join(final_metrics_path, metrics_txt), 'r')
    x = []
    y = []
    #full_xy = []
    for line in f:
        line = line.strip().split()
        x.append(float(line[0].strip())) 
        y.append(float(line[1].strip()))
    full_xy = [x, y]
    return full_xy

def scatter_plot(xy_list, txt_file):
    # visualize quality
    x = xy_list[0]
    y = xy_list[1]
    #print(y)
    plt.figure(figsize=(35, 10))
    plt.plot(x, y, 'o', color='black')
    y_label = ''
    if '__q' in txt_file:
        plt.ylim((0,93))
        y_label = 'Quality'
    elif '__sb' in txt_file:
        y_label = 'Strand Bias'
        plt.axhline(y=get_forward_strand_mean_std(xy_list), color='red')
        plt.axhline(y=get_reverse_strand_mean_std(xy_list), color='red')
    plt.title(y_label + ' at position x')
    plt.xlabel('Position')
    plt.ylabel(y_label)
    plt.savefig(os.path.basename(txt_file +'.png'))


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



# HEAT MAP BELOW #

# outputs [sample, gene, avgresidual]
def write_heat_map_avg_residual(txt_file, region_file, xy_list):
    # current sample
    sample_name = os.path.basename(txt_file)
    index_of_dot = sample_name.index('.')
    sample_name = sample_name[:index_of_dot]
    #print(sample_name)
    
    gene_name = os.path.basename(region_file)
    gene_name = os.path.splitext(gene_name)[0]
    #print(gene_name)    
    
    avg_res = get_avg_residual(xy_list)
    
    with open('heat_map_quality.txt', 'a') as Q_txt: 
        print(gene_name, '\t', sample_name, '\t', avg_res, '\t', file=Q_txt)
    Q_txt.close()
    
# outputs [sample, gene, proportion outside]
def write_heat_map_proportion_outside(txt_file, region_file, xy_list):
    # current sample
    sample_name = os.path.basename(txt_file)
    index_of_dot = sample_name.index('.')
    sample_name = sample_name[:index_of_dot]

    gene_name = os.path.basename(region_file)
    gene_name = os.path.splitext(gene_name)[0]
        
    proportion = get_proportion_outside(xy_list)
    
    with open('heat_map_strand_bias.txt', 'a') as SB_txt:
        print(gene_name, '\t', sample_name, '\t', proportion, '\t', file=SB_txt)
    SB_txt.close()

# heatmap metric for quality
def get_avg_residual(xy_list):
    y = xy_list[1]
    avg = np.average(y)
    all_residuals = []
    for i in y:
        residual = np.absolute(avg - i)
        all_residuals.append(residual)
    #print (np.average(all_residuals))
    return np.average(all_residuals)

# heatmap metric for strand bias
def get_proportion_outside(xy_list):
    # lower bound
    forward = get_forward_strand_mean_std(xy_list)
    # upper bound 
    reverse = get_reverse_strand_mean_std(xy_list)
    
    y = xy_list[1]
    
    outside = 0
    for i in y:
        if i < forward or i > reverse: outside += 1
    #print("strandbias: ", outside/len(y))
    try: proportion_outside = outside/len(y)
    except ZeroDivisionError: proportion_outside = -1
    return proportion_outside



def get_hm_data(heat_map_txt, num_samples):
    samples = []
    genes = []
    data = []
    #final = []
    with open(heat_map_txt, 'r') as hm_txt:
        count = 0
        row = []
        for l in hm_txt:
            A = l.rstrip().split()
            gene = A[0]
            sample = A[1]
            if gene not in genes: genes.append(gene)
            if sample not in samples: samples.append(sample)
            if count < num_samples: 
                row.append(float(A[2]))
                count += 1
            else:
                row = []
                row.append(float(A[2]))
                count = 1
            if count == num_samples-1: data.append(row)

    hm_txt.close()
    final = [genes, samples, data]
    return final
    
def plot_heat_map(num_samples):
    plot_path = '/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/chco-exome-analysis/hm_plots'
    #quality
    q_plot_name = 'quality_hm.png'
    q = get_hm_data('heat_map_quality.txt', num_samples)
    genes = q[0]
    samples = q[1]
    q_data = np.array(q[2])
    
    fig, ax = plt.subplots()
    #plt.figure(figsize=(30, 30))
    ax.set_xticks(np.arange(len(samples)))
    ax.set_yticks(np.arange(len(genes)))
    ax.set_xticklabels(samples)
    ax.set_yticklabels(genes)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.imshow(q_data)
    plt.colorbar(cmap='hot')
    
    # these for loops add values to the cells, but i would need to truncate and get them 
    #for i in range(len(genes)):
        #for j in range(len(samples)):
            #text = ax.text(j, i, q_data[i, j], ha="center", va="center", color="w")
    
    #plt.figure(figsize=(30,30))
    #ax.imshow(q_data, cmap='hot', interpolation='nearest')
    #plt.colorbar(cmap='hot')
    #ax.set_xticks(np.arange(len(samples)))
    #ax.set_yticks(np.arange(len(genes)))
    plt.savefig(os.path.join(plot_path, q_plot_name))
    #print (type(q_data), q_data)
    
    #strandBias
    sb = get_hm_data('heat_map_strand_bias.txt', num_samples)
    
    
    

#main()
