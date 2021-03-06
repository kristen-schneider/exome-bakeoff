import pysam
import os
from quality import main_quality


def main_single_pileup(tbx_pileup_file, regions_dir, metrics_list, sample_name):
    read_regions(tbx_pileup_file, regions_dir, metrics_list, sample_name)


def read_regions(tbx_pileup_file, regions_dir, metrics_list, sample_name):

    # regions_dir
    region = regions_dir.split('/')[7:9]

    print(region)

    # open files to store metrics for plotting
    for metric in metrics_list:
        metric_output_txt = open("/Users/krsc0813/exome-bakeoff/Analyses/" + metric + "/" + region[0] + "/" + region[1] + "/metric_files/" + sample_name + "_" + metric + ".txt", 'a')
        metric_output_txt.truncate(0)

    #regions_file_read = open(regions_dir, 'r')

    # iterate through all bed files in regions directory
    for bed_file in os.listdir(regions_dir):
        gene_bed = open(regions_dir + bed_file, 'r')
        gene_name = bed_file.split('.')[0]

        # go through the lines for each region file
        for line in gene_bed:
            line = line.rstrip().split()

            # chrm, start, end
            # handle X and Y chromosomes (strings)
            try: line_chrm = int(line[0])
            except ValueError: line_chrm = line[0]
            line_start = int(line[1])
            line_end = int(line[2])

            # find line (chrm, start, end) from gene.bed in sample
            for row in tbx_pileup_file.fetch(line_chrm, line_start, line_end, parser=pysam.asBed()):
                # get necessary values
                chrm = row[0]
                start = row[1]
                end = row[2]

                if "quality" in metrics_list:
                    quality_output_txt = open("/Users/krsc0813/exome-bakeoff/Analyses/" + metric + "/" + region[0] + "/" + region[1] + "/metric_files/" + sample_name + "_" + metric + ".txt", 'a')

                    quality = main_quality(row)
                    # format: chrm, start, end, quality_score
                    quality_output_txt.write(
                        chrm + '\t' + start + '\t' + end + '\t' + str(quality) + '\t' + gene_name + '\n')

                if "strandbias" in metrics_list:
                    strandbias = 0
                    # call strandbias script here

                if "noise" in metrics_list:
                    noise = 0
                    # call noise script here

                # ... and so on




    # for region in regions_file_read:
    #     region = region.rstrip().split()
    #     # convert chrm, start, end to proper type
    #     try: region_chrm = int(region[0])
    #     except ValueError: region_chrm = region[0]
    #     region_start = int(region[1])
    #     region_end = int(region[2])
    #
    #     # fetch regions
    #     for row in tbx_pileup_file.fetch(region_chrm, region_start, region_end, parser=pysam.asBed()):
    #         # get necessary values
    #         q_chrm = row[0]
    #         q_start = row[1]
    #         q_end = row[2]
    #         if "quality" in metrics_options:
    #             quality = main_quality(row)
    #             # format: chrm, start, end, quality_score
    #             quality_output_txt.write(q_chrm + '\t' + q_start + '\t' + q_end  + '\t' + str(quality) + '\t' + region[4] + '\n')
    #
    #         if "strandbias" in metrics_options:
    #             strandbias = 0
    #             # call strandbias script here
    #
    #         if "noise" in metrics_options:
    #             noise = 0
    #             # call noise script here
    #
    #         # ... and so on


