import numpy as np
import matplotlib.pyplot as plt
import os
import numpy as np

heatmap_metrics_path = '/Users/kristen/PycharmProjects/exome-bakeoff/txtFiles/heatmap_metrics/'

def boxplot_main():
    data=boxplot_data()
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    ax1.set_xticklabels(np.repeat(data[0], 2),
                        rotation=45, fontsize=8)
    
    #ax1.figure(figsize=(100, 30))
    ax1.boxplot(data[1])
    fig.savefig('boxplot.png', dpi=100)

def boxplot_data():
    sample_list = []
    all_qualities = []
    for sample in os.listdir(heatmap_metrics_path):
        sample_list.append(sample)
        curr_sample = open(heatmap_metrics_path+sample, 'r')
        curr_sample_range = []
        for line in curr_sample:
            line = line.rstrip().split()
            Q = float(line[1])
            curr_sample_range.append(Q)
        all_qualities.append(curr_sample_range)
    everything = [sample_list, all_qualities]
    return everything

boxplot_main()

# # Fixing random state for reproducibility
# np.random.seed(19680801)
#
# # fake up some data
# # Fixing random state for reproducibility
# np.random.seed(19680801)
#
# # fake up some data
# spread = np.random.rand(50) * 100
# center = np.ones(25) * 50
# flier_high = np.random.rand(10) * 100 + 100
# flier_low = np.random.rand(10) * -100
# data = np.concatenate((spread, center, flier_high, flier_low))
# spread = np.random.rand(50) * 100
# center = np.ones(25) * 40
# flier_high = np.random.rand(10) * 100 + 100
# flier_low = np.random.rand(10) * -100
# d2 = np.concatenate((spread, center, flier_high, flier_low))
# data.shape = (-1, 1)
# d2.shape = (-1, 1)
#
# data = [data, d2, d2[::2, 0]]
# fig7, ax7 = plt.subplots()
# ax7.set_title('Multiple Samples with Different sizes')
# ax7.boxplot(data)
# plt.savefig('boxplot.png', dpi=100)
