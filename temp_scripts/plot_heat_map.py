import sys
import os
import numpy as np
import matplotlib.pyplot as plt

data_files = sys.argv[1]
#quality_hm_data = sys.argv[1]
#strand_bias_hm_data = sys.argv[2]

num_samples = len([name for name in os.listdir('/scratch/Shares/layer/nextflow/kristen/fastq_to_vcf/mpileup/mpileup_output') if name.endswith('.tbi')])

def main():
    for data_file in os.listdir(data_files):
        if data_file.endswith('quality.txt'): plot_heat_map('Quality', data_file)
        if data_file.endswith('strand_bias.txt'): plot_heat_map('Strand Bias', data_file) 

def plot_heat_map(title, data_file):
    d = get_hm_data(data_file, num_samples)
    genes = d[0]
    samples = d[1]
    numeric_data = np.array(d[2])

    fig, ax = plt.subplots()
    fig.set_size_inches(35.0, 25.0)
    ax.set_xticks(np.arange(len(samples)))
    ax.set_yticks(np.arange(len(genes)))
    ax.set_xticklabels(samples)
    ax.set_yticklabels(genes)
    plt.title(title)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.imshow(numeric_data)
    plt.colorbar(cmap='hot')
    plt.savefig(title+'.png', dpi=100)

    # strand_bias 
    #sb_plot_name = 'strand_bias_hm.png'
    #sb = get_hm_data(strand_bias_hm_data, num_samples)
    #genes = sb[0]
    #samples = sb[1]
    #sb_data = np.array(q[2])

    #fig, ax = plt.subplots()
    #fig.set_size_inches(25.0, 15)
    #ax.set_xticks(np.arange(len(samples)))
    #ax.set_yticks(np.arange(len(genes)))
    #ax.set_xticklabels(samples)
    #ax.set_yticklabels(genes)
    #plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    #plt.imshow(sb_data)
    #plt.colorbar(cmap='hot')
    #plt.savefig(sb_plot_name, dpi=100)



def get_hm_data(heat_map_txt, num_samples):
    samples = []
    genes = []
    data = []
    #final = []
    with open(os.path.join(data_files, heat_map_txt), 'r') as hm_txt:
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
    
main()
