import argparse
import os
import pysam
from single_pileup import main_single_pileup

def main_pileup():
    # store arguments from command line
    args = get_cmdln_arguments()

    pileup_path = args['p']
    regions_path = args['r']
    metrics_options = args['m']

    # for each pileup file, call metrics script for all regions
    for pileup_file in os.listdir(pileup_path[0]):
        if (".bed" in pileup_file) and (".tbi" not in pileup_file):
            print(pileup_path[0] + pileup_file)
            sample_name = pileup_file.split('.')[0]
            # tabix the pileupfile for easy search
            tbx_pileup_file = pysam.TabixFile(pileup_path[0] + pileup_file)

            main_single_pileup(tbx_pileup_file, regions_path[0], metrics_options, sample_name)
            


# get the arguments from commandline run
def get_cmdln_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', nargs = 1, required = True, help = 'specify the directory of mpileup files')
    parser.add_argument('-r', nargs = 1, required = True, help = 'specify the directory of regions')
    parser.add_argument('-m', nargs = '+', required = True,  help = 'calculate specified metrics of samples')
    args = parser.parse_args()
    return vars(args)

main_pileup()
