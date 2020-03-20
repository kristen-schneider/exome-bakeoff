import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sns.set()
sns.set_style("white")

"""
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

# Ribbon plot of just the base quality
ax = sns.lineplot(x="GC", y="MEAN_BASE_QUALITY", data=data)
ax.legend_.remove()

# line plot of NORMALIZED_COVERAGE messy version
ax = sns.lineplot(x="GC", y="NORMALIZED_COVERAGE", hue="SAMPLE", data=data)
ax.legend_.remove()

# line plot of NORMALIZED_COVERAGE
ax = sns.lineplot(x="GC", y="NORMALIZED_COVERAGE", data=data)
ax.legend_.remove()

# bar plot
df['Normalize Windows'] = df['WINDOWS'] / sum(df['WINDOWS'])
ax = sns.barplot(x="GC", y="Normalize Windows", data=df, color='r').set(xticklabels=[x if x % 10 == 0 else None for x in range(1,100)])


# combined the plots
# rename column names
data.columns = ['ACCUMULATION_LEVEL', 'READS_USED', '% GC', 'WINDOWS', 'READ_STARTS',
                'Mean base quality', 'Normalized coverage', 'ERROR_BAR_WIDTH', 'SAMPLE',
                'LIBRARY', 'READ_GROUP']

# define figure size
fig = plt.figure(figsize=(10, 7))

# bar plot of windows, windows stacks are the same across all samples so only data from one
df['Normalize Windows'] = df['WINDOWS'] / sum(df['WINDOWS']) * 100
ax = sns.barplot(x="GC", y="Normalize Windows", data=df, color='r').set(xticklabels=[x if x % 10 == 0 else None for x in range(0,101)])

# make the two line plots
sns.lineplot(data=data, x='% GC', y='Normalized coverage', color="g")
ax2 = plt.twinx()
splot = sns.lineplot(data=data, x='% GC', y='Mean base quality', color="b", ax=ax2)
splot.set(ylim=(0, 30))

# define the colors and lines to be used in the legend
custom_lines = [Line2D([0], [0], color='g', lw=4),
                Line2D([0], [0], color='b', lw=4),
                Line2D([0], [0], color='r', lw=4)]

# create a custom legend
fig.legend(custom_lines, ['Normalized coverage', 'Mean base quality', 'Windows at % GC'], frameon=False)

# save the figure
plt.savefig('gc_bias_summary_plot.png')