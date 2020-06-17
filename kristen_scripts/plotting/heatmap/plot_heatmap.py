import matplotlib.pyplot as plt
import numpy as np
from read_data_sheet import read_data


def main():
    all_data = read_data('/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/plotting/csv_files/full_name_by_lp.csv')
    plot_heatmap(all_data)

    
def plot_heatmap(all_data):
    genes = all_data[0]
    sample_id = all_data[1]
    data = all_data[2]

    x_axis_label = genes
    y_axis_label = []
    for s in sample_id: y_axis_label.append('-'.join(s.split('-')[0:2]))

    plt.figure(figsize=(40, 50))
    plt.imshow(data,interpolation=None, cmap='coolwarm', vmin=35, vmax=58)
    plt.xticks(np.arange(0, len(genes)), x_axis_label, rotation=65)
    plt.xlabel('GENES')
    plt.yticks(np.arange(0, len(sample_id)), y_axis_label, rotation=0)
    plt.ylabel('SAMPLES')
    plt.colorbar()

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
    plt.savefig('full_name_by_lp.png')
    

main()
