import os

FULL_data='/Users/krsc0813/exome-bakeoff/kristen_scripts/full_name.csv'
LP_data='/Users/krsc0813/exome-bakeoff/kristen_scripts/lp_name.csv'
CAP_data='/Users/krsc0813/exome-bakeoff/kristen_scripts/capture_name.csv' 

#def read_data_sheet_main():
    
#    full_data = read_data(FULL_data)
#    lp_data = read_data(LP_data)
#    CAP_data = read_data(CAP_data)


def read_data(data_file):
    label_data = []
    metric_data = []
    genes = ''
    for line in open(data_file):
        A = line.split(',')
        
        # genes
        if len(genes) < 1: genes = A[1].rstrip().lstrip().split(' ')
        
        else:
            # name of full/lp/cap
            B = A[0].rstrip().lstrip()
            label_data.append(B)
            # data
            C = A[1].rstrip().lstrip().split(' ')
            data = [float(i) for i in C]
            metric_data.append(data)

    return ([genes, label_data, metric_data])
