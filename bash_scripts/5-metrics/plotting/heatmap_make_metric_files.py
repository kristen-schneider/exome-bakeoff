import os

samples_dir = '/Shares/layer_shared/projects/chco/exome_bakeoff/mpileup_downsample_tbi/'
regions_dir = '/Shares/layer_shared/projects/chco/exome_bakeoff/regions/59_ACMG/'
metrics_dir = '/scratch/Shares/layer/projects/chco/kristen/exome-bakeoff/Analysis/noise/59_ACMG/metric_files/'

def main():
    for m_file in os.listdir(metrics_dir):
        metrics_file = open(metrics_dir+m_file, 'r')
        sample_name = ('_').join(m_file.split('_')[0:3])
        heatmap_metrics_file = open('/scratch/Shares/layer/projects/chco/kristen/exome-bakeoff/Analysis/noise/59_ACMG/heatmap_files/'
            +sample_name+'_heatmap.txt', 'w')

        sample_dict = dict()
        header = ''
        for line in metrics_file:
            if header == '': header = line
            else:
                A = line.split()
                metric = A[3]
                gene = A[4]
                try: sample_dict[gene].append(float(metric))
                except KeyError: sample_dict[gene] = [float(metric)]
        #for d in sample_dict: print(d)
        for d in sample_dict.keys():
            #print(sum(sample_dict[d])/len(sample_dict[d]))
            heatmap_metrics_file.write(d+'\t'+str(sum(sample_dict[d])/len(sample_dict[d]))+'\n')
    
main()
