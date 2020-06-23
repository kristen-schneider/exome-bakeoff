import os
<<<<<<< HEAD
=======
from alphabetize import alphabetize_full
from alphabetize import alphabetize_prep
from alphabetize import alphabetize_capture
>>>>>>> 071cf845c105a150306b9a88322e497d5b1e93ba

heatmap_metrics_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/59_ACMG/downsampled/heatmap_metrics/'

def data_sheet_main():
    
    full_data = full_name_dict()
    lp_data = library_prep_dict(full_data)
    c_data = capture_dict(full_data)
    header = write_header()
    
    # lp-capture-sample
    full_name_data_file = open("full_name.csv", 'a')
    full_name_data_file.truncate(0)
    full_name_data_file.write(header + '\n')
    sorted_full = alphabetize_full(full_data)    

    for s in range(len(sorted_full)):
        full_name_data_file.write(str(sorted_full[s][0]) + ',' + (' ').join(map(str, sorted_full[s][1])).replace('[', '').replace(']', '') + '\n')
    
    # lp
    lp_data_file = open("lp_name.csv", 'a')
    lp_data_file.truncate(0)
    lp_data_file.write(header + '\n')
    sorted_prep = alphabetize_prep(lp_data)
    for p in range(len(sorted_prep)):
        for p2 in sorted_prep[p][1]:
            lp_data_file.write(str(sorted_prep[p][0]) + ',' + (' ').join(map(str,p2)).replace('[', '').replace(']', '') + '\n')
    #for a in sorted_prep:
    #    for b in sorted_prep[a]:
    #        lp_data_file.write(str(a) + ',' + (' ').join(map(str, b)) + '\n')
    #lp_data_file.write('\n\n')

    # capture
    capture_data_file = open("capture_name.csv", 'a')
    capture_data_file.truncate(0)
    capture_data_file.write(header + '\n')
    sorted_capture = alphabetize_capture(c_data)
    for c in range(len(sorted_capture)):
        for c2 in sorted_capture[c][1]:
            capture_data_file.write(str(sorted_capture[c][0]) + ',' + (' ').join(map(str,c2)).replace('[', '').replace(']', '') + '\n')
    #for a in sorted_capture:
    #    for b in sorted_capture[a]:
    #        capture_data_file.write(str(a) + ',' + (' ').join(map(str,b)) + '\n')
    #capture_data_file.write('\n\n')

def write_header():
    f = open('/Users/krsc0813/exome-bakeoff/Analyses/quality/59_ACMG/downsampled/heatmap_metrics/AgilentQXT-IDT-0720ME25_S12_L001_quality_hm.txt')
    
    genes_list = []
    genes_string = ''
    # genes as list
    for line in f:
        line = line.rstrip().split()
        genes_list.append(line[0])
    
    genes_string = (' ').join(genes_list)
    
    header = 'sample-name'+','+genes_string
    return header

def full_name_dict():
    sample_data = dict()
    for sample in os.listdir(heatmap_metrics_path):
        curr_sample = open(heatmap_metrics_path+sample, 'r')
        curr_sample_data = []
        for line in curr_sample:
            line = line.rstrip().split()
            Q = float(line[1])
            curr_sample_data.append(Q)
        sample_data[sample.split('_')[0]] = curr_sample_data

    return sample_data

def library_prep_dict(sample_data):
    library_prep_data = dict()
    for sample in sample_data:
        library_prep = sample.split('-')[0]+'-'+sample.split('-')[2]
        data = sample_data[sample]

        try: library_prep_data[library_prep].append(data)
        except KeyError: library_prep_data[library_prep] = [data]
    
    return library_prep_data

def capture_dict(sample_data):
    capture_data = dict()
    for sample in sample_data:
        capture = sample.split('-')[1]+'-'+sample.split('-')[2]
        data = sample_data[sample]

        try: capture_data[capture].append(data)
        except KeyError: capture_data[capture] = [data]

    return capture_data

        
data_sheet_main()
