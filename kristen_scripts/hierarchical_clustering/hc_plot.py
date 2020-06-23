import argparse
import pandas as pd
import numpy as np
import os
import seaborn as sns; sns.set(color_codes=True)
from matplotlib import pyplot as plt


def cluster_plot_main():
    args = get_cmdln_arguments()
    csv_file = args['csv'][0]
    cluster_plot(csv_file)
     
def cluster_plot(csv_file):
    csv = pd.read_csv(csv_file, index_col=0)

    # color by library prep
    library_prep = csv.pop('library prep tech')
    lut = dict(zip(np.unique(library_prep), 'mkygbr'))
    library_row_colors = pd.Series(library_prep, name='library prep tech').map(lut)
    
    # seaborn clustering
    g = sns.clustermap(csv,
                       vmin=35, vmax=58,
                       figsize=(25, 30),
                       cmap="coolwarm",
                       row_colors=library_row_colors)

    plt.savefig('59_clustering_lib_prep.png',
                dpi=150,
                figsize=(10, 10))

def get_cmdln_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-csv', nargs = 1, required = True, help = 'path to csv input file')
    args = parser.parse_args()
    return vars(args)

cluster_plot_main()

#import pandas as pd
#import matplotlib.patches as mpatches
#import matplotlib.pyplot as plt

#import numpy as np
#from matplotlib import pyplot as plt
#from hierarchical_make_csv_quality import sample_gene_quality_technology
#from hierarchical_make_csv_quality import gene_technology

#from scipy.cluster.hierarchy import dendrogram
#from sklearn.datasets import load_iris
#from sklearn.cluster import AgglomerativeClustering


## https://cmdlinetips.com/2020/01/heatmaps-with-seaborns-clustermap/

