import os
import sys
import pysam


sample_name = sys.argv[1]
pileup_file = sys.argv[2]
VCF_dir = sys.argv[3]
regions = sys.argv[4]
out_dir = sys.argv[5]

def test():
    # using pysam to filter through pileup
    pileup_tbx = pysam.TabixFile(pileup_file)
    for vcf_file in os.listdir(VCF_dir):
        if sample_name in vcf_file and ".gz" in vcf_file and ".tbi" not in vcf_file:
            vcf_tbx = pysam.TabixFile(VCF_dir+vcf_file)
    
    
    #test_gene = "LDLR.bed"
    #LDLR_bed = open(regions+test_gene, 'r')
    for gene_bed in os.listdir(regions):
        gene_name = gene_bed.split('.')[0]
        
        for line in open(regions+gene_bed, 'r'):
            #print(line)
            line = line.rstrip().split()
            try: line_chrm = int(line[0])
            except ValueError: line_chrm = line[0]
            line_start = int(line[1])
            line_end = int(line[2])
    
            for pileup_hit in pileup_tbx.fetch(line_chrm, line_start, line_end, parser=pysam.asBed()):
                #print(pileup_hit)
                pileup_line = str(pileup_hit).split()
                try: pileup_chrm = int(pileup_line[0])
                except ValueError: pileup_chrm = pileup_line[0]
                pileup_start = int(pileup_line[1])
                pileup_end = int(pileup_line[2])
                pileup_read_count = pileup_line[4] 
                pileup_reads = pileup_line[5]
                print('pileup' + '\t'+ str(pileup_chrm) + '\t' + str(pileup_start) + '\t' + str(pileup_end) + '\t' + pileup_read_count + '\t' + pileup_reads + '\t' + gene_name)
                
                # fetch vcf  
                for vcf_hit in vcf_tbx.fetch(pileup_chrm, pileup_start, pileup_end, parser=pysam.asBed()):
                    vcf_line = str(vcf_hit).split()
                    try: vcf_chrm = int(vcf_line[0])
                    except ValueError: vcf_chrm = vcf_line[0]
                    vcf_start = int(vcf_line[1])
                    vcf_end = int(vcf_start + 1)
                    vcf_ref = vcf_line[3]
                    vcf_alt = vcf_line[4]
                    vcf_genotype = vcf_line[9].split(':')[0]
                    #print('vcf' + '\t' + str(vcf_line))
                    print('vcf' + '\t'+ str(vcf_chrm) + '\t' + str(vcf_start) + '\t' + str(vcf_end) + '\t' + vcf_ref + '\t' + vcf_alt + '\t' + vcf_genotype + '\t' + gene_name)

    pileup_tbx.close()
    #regions+gene_bed.close()        
  
                
def main():
    # using pysam to filter through pileup
    pileup_tbx = pysam.TabixFile(pileup_file)
    
    # iterate through all gene.bed files in regions dir
    for gene_bed in os.listdir(regions):
        # gene name
        gene_name = gene_bed.split('.')[0]
        # open and read bed file
        gene_bed_open = open(regions+gene_bed, 'r')
        # get info for regions
        print(gene_name)
        for line in gene_bed_open:
            line = line.rstrip().split()
            try: line_chrm = int(line[0])
            except ValueError: line_chrm = line[0]
            line_start = int(line[1])
            line_end = int(line[2])
            
            
            #found_in_VCF = False

            # pass in regions to vcf files
            for vcf_file in os.listdir(VCF_dir):
                print(vcf_file)
                # make sure to only target the proper vcf file
                if ".gz" in vcf_file and ".tbi" not in vcf_file:
                    vcf_tbx = pysam.TabixFile(VCF_dir+vcf_file)
                    for vcf_hit in vcf_tbx.fetch(line_chrm, line_start, line_end, parser=pysam.asBed()):
                        # store interesting values from vcf file
                        try: vcf_chrm = int(vcf_hit[0])
                        except ValueError: vcf_chrm = vcf_hit[0]
                        vcf_pos = int(vcf_hit[1])
                        vcf_ref = vcf_hit[3]
                        vcf_alt = vcf_hit[4]
                        vcf_genotype = vcf_hit[9].split(':')[0]
                        
                        print(vcf_genotype)
                        #found_in_VCF = True
                    #if found_in_VCF:
                        #print("alt")
                        #found_in_VCF = False
                    #else:
                        #print("homo_ref")

test()

