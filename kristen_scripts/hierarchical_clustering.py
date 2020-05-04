import os
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering


## https://cmdlinetips.com/2020/01/heatmaps-with-seaborns-clustermap/

TXTFILES = '/Users/kristen/PycharmProjects/exome-bakeoff/txtFiles/heatmap_metrics/'
panda_file = '/Users/kristen/PycharmProjects/exome-bakeoff/quality-pandas.csv'
test_array = [[0,0,0],[0,0,0],[0,0,0],[1,1,1],[0,1,1],[0,0,1]]


## SEABORN
def seaborn_main():
    make_pandas_file()
    gapminder = read_pandas()


    #twoD_array = make_2D_array()
    #seaborn_plot(twoD_array)

def read_pandas():
    data_url = 'http://bit.ly/2cLzoxH'
    # read data from url as pandas dataframe
    gapminder = pd.read_csv(panda_file)
    print(gapminder.head(3))
    heatmap_data = pd.pivot_table(gapminder, values='quality',
                                  index=['sample'],
                                  columns='gene')
    print(heatmap_data.head())
    sns.clustermap(heatmap_data)
    plt.savefig('heatmap_with_Seaborn_clustermap_python.png',
                dpi=150, figsize=(8, 12))

def make_pandas_file():
    pandas_file = open('quality-pandas.csv', 'a')
    pandas_file.truncate(0)
    pandas_file.write('sample' + ',' + 'gene' + ',' + 'quality' + '\n')

    for sample in os.listdir(TXTFILES):
        sample_name = ('-'.join(sample.split('-')[0:2]))

        s = open(TXTFILES + sample, 'r')
        for line in s:
            A = line.rstrip().split()
            gene = A[0]
            quality = A[1]
            pandas_file.write(sample_name + ',' + gene + ',' + quality + '\n')


seaborn_main()

# OLD CODE FOR 2D ARRAY CONSTRUCTION
# This way does not include any labeling
#
# def seaborn_plot(twoD_array):
#     #iris = sns.load_dataset("iris")
#     #species = iris.pop("species")
#     g = sns.clustermap(test_array)
#     g.savefig("dendrogram")
#
# def make_2D_array():
#     twoD_array = []
#
#     for sample in os.listdir(TXTFILES):
#         sample_array = make_sample_array(sample)
#         twoD_array.append(sample_array)
#     return twoD_array
#
# def make_sample_array(sample):
#     sample_array = []
#     sample_file = open(TXTFILES + sample, 'r')
#     for line in sample_file:
#         line = line.rstrip().split()
#         quality_score = float(line[1])
#         sample_array.append(quality_score)
#     return sample_array
#
#
#
