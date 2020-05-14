# This script will take one input file (bam or pileup)...
#...and pull it through all regions files
#...and compute whatever metrics are requested of it
# INPUT:(1) single bam OR single pileup
#       (2) directory of regions files

import pysam
from quality import main_quality


def main_single_pileup(tbx_pileup_file, regions_file, metrics_options, sample_name):
    read_regions(tbx_pileup_file, regions_file, metrics_options, sample_name)


def read_regions(tbx_pileup_file, regions_file, metrics_options, sample_name):

    # open files to store metrics for plotting
    quality_output_txt = open("/Users/kristen/PycharmProjects/exome-bakeoff/Analyses/quality/metric_files/" + sample_name + "-QUALITY.txt", 'a')
    quality_output_txt.truncate(0)

    regions_file_read = open(regions_file, 'r')

    for region in regions_file_read:
        region = region.rstrip().split()
        # convert chrm, start, end to proper type
        try: region_chrm = int(region[0])
        except ValueError: region_chrm = region[0]
        region_start = int(region[1])
        region_end = int(region[2])

        # fetch regions
        for row in tbx_pileup_file.fetch(region_chrm, region_start, region_end, parser=pysam.asBed()):
            # get necessary values
            q_chrm = row[0]
            q_start = row[1]
            q_end = row[2]
            if "quality" in metrics_options:
                quality = main_quality(row)
                # format: chrm, start, end, quality_score
                quality_output_txt.write(q_chrm + '\t' + q_start + '\t' + q_end  + '\t' + str(quality) + '\t' + region[4] + '\n')

            if "strandbias" in metrics_options:
                strandbias = 0
                # call strandbias script here

            if "noise" in metrics_options:
                noise = 0
                # call noise script here

            # ... and so on


