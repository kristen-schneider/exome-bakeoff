import matplotlib.pyplot as plt
import numpy as np
from heatmap_read_csv import read_csv_main


def main():
    all_data = read_csv_main()
    plot_heatmap(all_data)
    #all_data = read_data('/Users/krsc0813/exome-bakeoff/kristen_scripts/plotting/csv_files/59_ACMG/full_name_by_lp.csv')
    #plot_heatmap(all_data)


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
    plt.imshow(data, interpolation=None, cmap='Blues', vmin=0, vmax=1)
    plt.xticks(np.arange(0, len(genes)), x_axis_label, rotation=65)
    plt.xlabel('GENES')
    plt.yticks(np.arange(0, len(sample_id)), y_axis_label, rotation=0)
    plt.ylabel('SAMPLES')
    plt.colorbar()
    
    plt.savefig('noise_59ACMG.png')
    
    # plt.setp(rotation=45, ha="right", rotation_mode="anchor")
    # plt.colorbar()
    # plt.plot([], c='green', label='Agilent')
    # plt.plot([], c='blue', label='IDT')
    # plt.plot([], c='red', label='Roche')
    # plt.plot([], c='yellow', label='Twist')
    # plt.legend()
    #
    # plt.title('Quality by Library Prep\n(colored by Capture Technology)')
    # plt.xticks(np.arange(1, len(sample_data[1]) + 1), x_label, rotation=65)
    # plt.xlim(-1, len(sample_data[1]) + 1)

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
