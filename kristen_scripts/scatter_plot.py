sample_qualities = '/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/metric_files'
for sample in sample_qualities:
    curr_sample = open(sample, 'r')
    for line in curr_sample:
        continue