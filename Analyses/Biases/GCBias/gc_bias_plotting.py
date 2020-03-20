import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re

sns.set()
sns.set_style("white")

"""
import os

os.chdir('Analyses/Biases/GCBias')
"""

data = None
first = True
first_df = None
for d in os.listdir('results'):
    for f in os.listdir('results/' + d):
        if f == 'gc_bias_metrics.txt':
            df = pd.read_csv('results/' + d + '/gc_bias_metrics.txt', sep='\t', skiprows=6)
            df['SAMPLE'] = d
            if first:
                data = df
                first = False
                first_df = df
            else:
                data = data.append(df)

# make a messy explority plot of all the base qualities plotted individually
ax = sns.lineplot(x="GC", y="MEAN_BASE_QUALITY", hue="SAMPLE", data=data)
ax.legend_.remove()
plt.clf()

# Ribbon plot of just the base quality
ax = sns.lineplot(x="GC", y="MEAN_BASE_QUALITY", data=data)
plt.clf()

# line plot of NORMALIZED_COVERAGE messy version
ax = sns.lineplot(x="GC", y="NORMALIZED_COVERAGE", hue="SAMPLE", data=data, color='b')
ax.legend_.remove()
plt.savefig('normalized_coverage_messy.png')
plt.clf()

# line plot of NORMALIZED_COVERAGE
ax = sns.lineplot(x="GC", y="NORMALIZED_COVERAGE", data=data)
plt.clf()

# bar plot
df['Normalize Windows'] = df['WINDOWS'] / sum(df['WINDOWS'])
ax = sns.barplot(x="GC", y="Normalize Windows", data=df, color='r').set(
    xticklabels=[x if x % 10 == 0 else None for x in range(1, 100)])
plt.clf()

# combined the plots
# rename column names
data.columns = ['ACCUMULATION_LEVEL', 'READS_USED', '% GC', 'WINDOWS', 'READ_STARTS',
                'Mean base quality', 'Normalized coverage', 'ERROR_BAR_WIDTH', 'SAMPLE',
                'LIBRARY', 'READ_GROUP']

# define figure size
fig = plt.figure(figsize=(10, 7))

# bar plot of windows, windows stacks are the same across all samples so only data from one
df['Normalize Windows'] = df['WINDOWS'] / sum(df['WINDOWS']) * 100
ax = sns.barplot(x="GC", y="Normalize Windows", data=df, color='r').set(
    xticklabels=[x if x % 10 == 0 else None for x in range(0, 101)])

# make the two line plots
sns.lineplot(data=data, x='% GC', y='Normalized coverage', color="b", ci='sd')
ax2 = plt.twinx()
splot = sns.lineplot(data=data, x='% GC', y='Mean base quality', color="g", ax=ax2, ci='sd')
splot.set(ylim=(0, 35))

# define the colors and lines to be used in the legend
custom_lines = [Line2D([0], [0], color='g', lw=4),
                Line2D([0], [0], color='b', lw=4),
                Line2D([0], [0], color='r', lw=4)]

# create a custom legend
fig.legend(custom_lines, ['Mean base quality', 'Normalized coverage', 'Windows at % GC'], frameon=False)

# save the figure
plt.savefig('gc_bias_summary_plot.png')
plt.clf()

def make_plot_group_by_name_2():
    plt.figure(figsize=(10, 7))
    data = None
    first = True
    groups = []
    for d in os.listdir('results'):
        groups.append(re.sub('-.*', '', re.sub('^\\w+-', '', d)))
    colors = ['r', 'b', 'g', 'purple', 'orange']
    groups = list(set(groups))
    color_map = { groups[i]:colors[i] for i in range(len(groups))}
    print(color_map)
    for d in os.listdir('results'):
        group = re.sub('-.*', '', re.sub('^\\w+-', '', d))
        for f in os.listdir('results/' + d):
            if f == 'gc_bias_metrics.txt':
                df = pd.read_csv('results/' + d + '/gc_bias_metrics.txt', sep='\t', skiprows=6)
                plt.plot("GC", 'NORMALIZED_COVERAGE', data=df, color=color_map[group])
                # plt.plot("GC", "MEAN_BASE_QUALITY", data=df, color='g')
                df['SAMPLE'] = d
                if first:
                    data = df
                    first = False
                else:
                    data = data.append(df)

    # define the colors and lines to be used in the legend
    custom_lines = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4),
                    Line2D([0], [0], color=colors[3], lw=4)]

    # create a custom legend
    plt.legend(custom_lines, groups, frameon=False)
    plt.xlabel('% GC')
    plt.ylabel('Normalized coverage')
    plt.savefig('normalized_coverage_tech_2.png')
    plt.clf()

make_plot_group_by_name_2()


def make_plot_group_by_name_1():
    plt.figure(figsize=(10, 7))
    data = None
    first = True
    groups = []
    for d in os.listdir('results'):
        print(re.sub('-.*-.*', '', d))
        groups.append(re.sub('-.*-.*', '', d))
    colors = ['r', 'b', 'g', 'purple', 'orange','y']
    groups = list(set(groups))
    print(groups)
    color_map = { groups[i]:colors[i] for i in range(len(groups))}
    print(color_map)

    for d in os.listdir('results'):
        group = re.sub('-.*-.*', '', d)
        for f in os.listdir('results/' + d):
            if f == 'gc_bias_metrics.txt':
                df = pd.read_csv('results/' + d + '/gc_bias_metrics.txt', sep='\t', skiprows=6)
                plt.plot("GC", 'NORMALIZED_COVERAGE', data=df, color=color_map[group])
                # plt.plot("GC", "MEAN_BASE_QUALITY", data=df, color='g')
                df['SAMPLE'] = d
                if first:
                    data = df
                    first = False
                else:
                    data = data.append(df)

    # define the colors and lines to be used in the legend
    custom_lines = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4),
                    Line2D([0], [0], color=colors[3], lw=4),
                    Line2D([0], [0], color=colors[4], lw=4),
                    Line2D([0], [0], color=colors[5], lw=4)]

    # create a custom legend
    plt.legend(custom_lines, groups, frameon=False)
    plt.xlabel('% GC')
    plt.ylabel('Normalized coverage')
    plt.savefig('normalized_coverage_tech_1.png')
    plt.clf()

make_plot_group_by_name_1()