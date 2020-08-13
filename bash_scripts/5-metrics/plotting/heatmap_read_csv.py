import os

in_file = '/Users/krsc0813/exome-bakeoff/bash_scripts/5-metrics/plotting/noise_59ACMG.csv'

def read_csv_main():
    
    csv_f = open(in_file, 'r')
    genes_list = ''
    samples_list = []
    data = []
    
    for line in csv_f:
        if genes_list == '': genes_list = line.replace('sample name,', '').split(',')
        else:
            A = line.rstrip().split(',')
            samples_list.append(A[0])
            sample_data = A[1:]
            data.append(sample_data)    
             
    return [genes_list, samples_list, data]

#read_csv_main()
