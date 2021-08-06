import os

heatmap_metric_dir='/scratch/Shares/layer/projects/chco/kristen/exome-bakeoff/Analysis/noise/59_ACMG/heatmap_files/'
samples_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_downsample_tbi/'
genes_dir='/Shares/layer_shared/projects/chco/exome_bakeoff/regions/59_ACMG/'

def main():
    samples_list = get_samples_list(samples_dir)
    genes_list = get_genes_list(genes_dir)
    
    #key = sample name
    #value = [metric values in list in order of genes]
    heatmap_dict = dict()
        
    for txt in os.listdir(heatmap_metric_dir):
        dict_key = ('_').join(txt.split('_')[0:3])
        dict_value = get_sample_scores(txt, genes_list)
        heatmap_dict[dict_key] = dict_value
    
    # write csv
    write_csv(samples_list, genes_list, heatmap_dict)
        #sample_txt = open(heatmap_metric_dir+txt, 'r')
        #for line in sample_txt:


    #return [samples_list, genes_list, heatmap_dict

def write_csv(samples_list, genes_list, heatmap_dict):
    # TO CHANGE
    csv_f = open('noise_59ACMG.csv', 'w')
    
    # header
    csv_f.write('sample name')
    for g in genes_list: csv_f.write(',' + g)
    csv_f.write('\n')
    
    # data
    for s in heatmap_dict.keys():
        csv_f.write(s)
        str_values = str(heatmap_dict[s]).replace('[', '').replace(']','')
        csv_f.write(','+str_values+'\n')

# returns a list of all the scores for one sample in the order of genes listed
def get_sample_scores(txt, genes_list):
    sample_scores = ['NaN'] * len(genes_list)
    
    #sample_txt = open(heatmap_metric_dir+txt, 'r')
    #count = 0
    #for line in sample_txt: count+=1
    #sample_txt.close()

    sample_txt = open(heatmap_metric_dir+txt, 'r')
    for line in sample_txt:
        A = line.rstrip().split()
        gene = A[0]
        score = A[1]
            
        if gene in genes_list:
            index = genes_list.index(gene)
            sample_scores[index] = float(score)

    return sample_scores
  

def get_samples_list(samples_dir):
    samples_list = []
    for tbi in os.listdir(samples_dir):
        if ".txt" in tbi:
            sample = tbi.split('.')[0]
            samples_list.append(sample)
        else: continue
    return samples_list

def get_genes_list(genes_dir):
    genes_list = []
    for g in os.listdir(genes_dir):
        if ".bed" in g:
            gene = g.split('.')[0]
            genes_list.append(gene)
        else: continue
    return genes_list

main()
