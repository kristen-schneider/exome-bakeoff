import os

TXTFILES = '/Users/kristen/PycharmProjects/exome-bakeoff/txtFiles/heatmap_metrics/'
sample_gene_quality_technology = '/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/sample_gene_quality_technology-onco.csv'
gene_technology = '/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/gene_technology-onco.csv'
gene_sample = '/Users/kristen/PycharmProjects/exome-bakeoff/kristen_scripts/gene_sample-onco.csv'
Onco_genes_txt_file = '/Users/kristen/PycharmProjects/exome-bakeoff/regions/oncogenes.txt'
ACMG_genes_txt_file = '/Users/kristen/PycharmProjects/exome-bakeoff/regions/59_acmg_genes.txt'

def ACMG_genes():
    gene_names = []
    ACMG_genes_txt = open(ACMG_genes_txt_file, 'r')
    for line in ACMG_genes_txt:
        if 'SDHD' not in line: gene_names.append(line.rstrip())
    return gene_names


def sample_gene_quality_technology():
    sample_gene_quality_technology_csv = open('sample_gene_quality_technology.csv', 'a')
    sample_gene_quality_technology_csv.truncate(0)
    sample_gene_quality_technology_csv.write('sample' + ',' + 'gene' + ',' + 'quality' + ',' + 'technology' + '\n')

    for sample in os.listdir(TXTFILES):
        sample_name = ('-'.join(sample.split('-')[0:3]))

        technology = sample.split('-')[1]
        #technology = ('-'.join(sample.split('-')[0:2]))

        s = open(TXTFILES + sample, 'r')
        for line in s:
            A = line.rstrip().split()
            gene = A[0]
            quality = A[1]
            sample_gene_quality_technology_csv.write(sample_name + ',' + gene + ',' + quality + ',' + technology + '\n')


def gene_technology():
    # header
    gene_names = ACMG_genes()
    gene_technology_csv = open('gene_sample.csv', 'a')
    gene_technology_csv.truncate(0)
    for gene in gene_names: gene_technology_csv.write(gene + ',')
    gene_technology_csv.write('capture' + '\n')

    # body
    for sample in os.listdir(TXTFILES):
        # make heatmap_metric file into dictionary for easy search by gene
        sample_dict = heatmap_metric_to_dict(sample)
        for gene in gene_names:
            try: gene_technology_csv.write(sample_dict[gene] + ',')
            except KeyError: continue

        # technology name
        technology_name = sample.split('_')[0].split('-')[0]
        #technology_name = "-".join(technology_name)
        gene_technology_csv.write(technology_name + '\n')

def heatmap_metric_to_dict(sample):
    sample_dict = {}
    sample_txt = open(TXTFILES+sample, 'r')
    for line in sample_txt: sample_dict[line.rstrip().split()[0]] = line.rstrip().split()[1]
    return sample_dict


gene_technology()
sample_gene_quality_technology()
