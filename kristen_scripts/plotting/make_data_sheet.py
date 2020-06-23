import os

heatmap_metrics_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/59_ACMG/downsampled/heatmap_metrics/'



def data_sheet_main():
    
    full_data = full_name_dict()
    lp_data = library_prep_dict(full_data)
    c_data = capture_dict(full_data)
    
    # lp-capture-sample
    full_name_data_file = open("full_name.txt", 'a')
    full_name_data_file.truncate(0)
    for s in full_data: full_name_data_file.write(str(s) + ': ' + str(full_data[s]) + '\n')
    
    # lp
    lp_data_file = open("lp_name.txt", 'a')
    lp_data_file.truncate(0)
    for a in lp_data:
        for b in lp_data[a]:
            lp_data_file.write(str(a) + ': ' + str(b) + '\n')
    lp_data_file.write('\n\n')

    # capture
    capture_data_file = open("capture_name.txt", 'a')
    capture_data_file.truncate(0)
    for a in c_data:
        for b in c_data[a]:
            capture_data_file.write(str(a) + ': ' + str(b) + '\n')
    capture_data_file.write('\n\n')
    

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
        library_prep = sample.split('-')[0]
        data = sample_data[sample]

        try: library_prep_data[library_prep].append(data)
        except KeyError: library_prep_data[library_prep] = [data]
    
    return library_prep_data

def capture_dict(sample_data):
    capture_data = dict()
    for sample in sample_data:
        capture = sample.split('-')[1]
        data = sample_data[sample]

        try: capture_data[capture].append(data)
        except KeyError: capture_data[capture] = [data]

    return capture_data

        
data_sheet_main()
