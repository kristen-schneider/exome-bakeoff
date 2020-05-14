import os
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import numpy as np
from matplotlib import pyplot as plt
from hierarchical_make_csv_quality import sample_gene_quality_technology
from hierarchical_make_csv_quality import gene_technology

from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering


## https://cmdlinetips.com/2020/01/heatmaps-with-seaborns-clustermap/

TXTFILES = '/Users/kristen/PycharmProjects/exome-bakeoff/txt/heatmap_metrics/'
sample_gene_quality_technology = '/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/sample_gene_quality_technology.csv'
gene_technology = '/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/gene_technolgoy.csv'
gene_sample = '/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/gene_sample.csv'
test_array = [[0,0,0],[0,0,0],[0,0,0],[1,1,1],[0,1,1],[0,0,1]]


## SEABORN
def seaborn_main():
    #with_sample_labels()
    with_technolgoy_colors()

    # iris = sns.load_dataset("iris")
    # print(iris.head(50))

    #gapminder = read_pandas()

def with_sample_labels():
    csv = pd.read_csv(sample_gene_quality_technology)
    print(csv.head(3))

    technology = csv['technology']

    pal = sns.light_palette('green', np.unique(technology).size)
    lut = dict(zip(np.unique(technology), pal))
    row_colors = csv[technology].map(lut)


    heatmap_data = pd.pivot_table(csv,
                                  values='quality',
                                  index=['sample'],
                                  columns='gene')

    g = sns.clustermap(heatmap_data,
                       figsize=(25, 30),
                       cmap="RdBu",
                       row_colors=row_colors)
    g.fig.suptitle('QUALITY FOR DIFFERENT TECHNOLOGIES')

    plt.savefig('sample_labels.png',
                dpi=150, figsize=(10, 10))


def with_technolgoy_colors():
    csv = pd.read_csv(gene_sample)
    technology = csv.pop('capture')

    pal = sns.light_palette('green', np.unique(technology).size)
    lut = dict(zip(np.unique(technology), 'mkygbr'))
    row_colors = pd.Series(technology, name='library prep tech').map(lut)

    g = sns.clustermap(csv,
                       figsize=(25, 30),
                       cmap="RdBu",
                       row_colors=row_colors)
                       #col_cluster=False)
    g.fig.suptitle('QUALITY FOR SAMPLES\n(green column separation by six capture)')

    #g.set(xlabel='my x label', ylabel='my y label')


    plt.savefig('sample_colors.png',
                dpi=150,
                figsize=(10, 10))

seaborn_main()