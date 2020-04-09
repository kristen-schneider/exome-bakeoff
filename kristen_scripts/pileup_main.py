# This script will call "padre" (to run on each pileup file)
# INPUT:(1) directory of pileup files
#       (2) request for which metric we need
# --> for each pileup:
#   --> call "padre" script
#   --> "padre" script will call the metric script and write some metric files

import argparse
import os
import pysam
from single_pileup import main_single_pileup

def main_pileup():
    # store arguments from command line
    args = get_cmdln_arguments()

    pileup_path = args['p']
    regions_file = args['r']
    metrics_options = args['m']

    # for each pileup file, call metrics script for all regions
    for pileup_file in os.listdir(pileup_path[0]):
        if (".bed" in pileup_file) and (".tbi" not in pileup_file):

            # tabix the pileupfile for easy search
            tbx_pileup_file = pysam.TabixFile(pileup_path[0] + '/' + pileup_file)

            main_single_pileup(tbx_pileup_file, regions_file[0], metrics_options)
            


# get the arguments from commandline run
def get_cmdln_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', nargs = 1, required = True, help = 'specify the directory of mpileup files')
    parser.add_argument('-r', nargs = 1, required = True, help = 'specify the directory of regions')
    parser.add_argument('-m', nargs = '+', required = True,  help = 'calculate specified metrics of samples')
    args = parser.parse_args()
    return vars(args)

main_pileup()
