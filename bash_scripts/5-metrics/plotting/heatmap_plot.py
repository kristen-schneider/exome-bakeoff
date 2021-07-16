import matplotlib.pyplot as plt
import numpy as np
from heatmap_read_csv import read_csv_main


def main():
    all_data = read_csv_main()
    plot_heatmap(all_data)


def plot_heatmap(all_data):
    genes = all_data[0]
    sample_id = all_data[1]
    data = convert_data_to_float(all_data[2])

    #print(genes)
    #print(sample_id)
    #print(data)
    x_axis_label = genes
    y_axis_label = []
    for s in sample_id: y_axis_label.append('-'.join(s.split('-')[0:2]))
    

    plt.figure(figsize=(70, 30))
    plt.imshow(data, interpolation=None, cmap='Blues', vmin=-0.005, vmax=0)
    plt.xticks(np.arange(0, len(genes)), x_axis_label, rotation=65)
    plt.xlabel('GENES')
    plt.yticks(np.arange(0, len(sample_id)), y_axis_label, rotation=0)
    plt.ylabel('SAMPLES')
    plt.colorbar()
    
    plt.savefig('noise_59ACMG.png')
    

def convert_data_to_float(data):
    float_data = []
    for sample in data:
        float_sample = []
        for d in sample: 
            try: float_sample.append(float(d))
            except ValueError: float_sample.append(np.nan)#print(d)#continue#t_sample.append(d)
        float_data.append(float_sample)

    return float_data


main()
