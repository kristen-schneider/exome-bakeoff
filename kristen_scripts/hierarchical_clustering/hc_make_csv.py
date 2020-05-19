import os
import argparse

def csv_main():
    # store arguments from command line
    args = get_cmdln_arguments()
    heatmap_files = args['hm'][0]
    regions_file = args['r'][0]
    csv_path = args['o'][0]
    title = args['t'][0]
    
    write_csv(heatmap_files, regions_file, csv_path, title)

# return a list of gene names for lableing purposes (header for rows)
def get_gene_names(regions_file):
    gene_names = []
    gene_regions_txt = open(regions_file, 'r')
    for line in gene_regions_txt:
        # SDHD gene not found in Homo_sapiens.GRCh37.82.exons.bed
        if 'SDHD' not in line: gene_names.append(line.rstrip())
    return gene_names

# columns are genes
# rows are sample
# last column is a characteristic upon which we want to visualize clustering (e.g. library prep)
def write_csv(heatmap_files, regions_file, csv_path, title):
    # header
    gene_names = get_gene_names(regions_file)
    csv_file = open(csv_path + title + '.csv', 'a')
    csv_file.truncate(0)
    for gene in gene_names: csv_file.write(gene + ',')
    csv_file.write('library prep tech' + '\n')
    print(csv_file)
    # body
    for sample in os.listdir(heatmap_files):
        # make heatmap_metric file into dictionary for easy search by gene
        curr_sample_dict = heatmap_metric_to_dict(sample, heatmap_files)
        for gene in gene_names:
            
            try: csv_file.write(curr_sample_dict[gene] + ',')
            except KeyError: print(gene)

        # technology name
        library_prep_name = sample.split('_')[0].split('-')[0]
        #technology_name = "-".join(technology_name)
        csv_file.write(library_prep_name + '\n')

# make heatmap_metric file into dictionary for easy search by gene
def heatmap_metric_to_dict(sample, heatmap_files):
    sample_dict = {}
    sample_txt = open(heatmap_files + sample, 'r')
    for line in sample_txt: sample_dict[line.rstrip().split()[0]] = line.rstrip().split()[1]
    return sample_dict

# get the arguments from commandline run
def get_cmdln_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hm', nargs = 1, required = True, help = 'give location for heatmap metric files')
    parser.add_argument('-r', nargs = 1, required = True, help = 'give location for the regions text file (list)')
    parser.add_argument('-o', nargs = 1, required = True, help = 'give the output location for the csv file')
    parser.add_argument('-t', nargs = 1, required = True, help = 'give the plot a title')
    args = parser.parse_args()
    return vars(args)

csv_main()
