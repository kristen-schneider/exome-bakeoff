import os
import sys
import pysam


pileup_file = sys.argv[1]
VCF_dir = sys.argv[2]
region_dir = sys.argv[3]
#out_dir = sys.argv[4]

def test():
    for i in range(len(sys.argv)):
        print(sys.argv[i])

def main():
    print(pileup_file)
    # cycle through each gene in region file
    for region_file in os.listdir(region_dir):
        region_bed = open(region_dir+region_file, 'r')
        for line in region_bed:
            line = line.rstrip().split()
            try: line_chrm = int(line[0])
            except ValueError: line_chrm = line[0]
            line_start = int(line[1])
            line_end = int(line[2])        
    
        # cycle through each VCF per gene
        for vcf_file in os.listdir(VCF_dir):
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
                    # convert '0/1' to ['A','G']
                    allele_genotype = get_allele_genotype(vcf_genotype, vcf_ref, vcf_alt)
                
                    # find matching hits in pileup file
                    pileup_tbx = pysam.TabixFile(pileup_file)
                    for pileup_hit in pileup_tbx.fetch(vcf_chrm, vcf_pos, vcf_pos+1, parser=pysam.asBed()):
                        pileup_ref = pileup_hit[3]
                        pileup_read_count = float(pileup_hit[4])
                        pileup_base_read = pileup_hit[5]
                        ref_matches = find_ref_matches(pileup_base_read)
                        alt_matches = find_alt_matches(pileup_base_read, allele_genotype)

                        try: match_contribution = float(ref_matches+alt_matches)/float(pileup_read_count)
                        except ZeroDivisionError: continue
                        print(str(vcf_chrm) + '\t' + str(vcf_pos) + '\t' + str(vcf_pos+1) + 
                            '\t' + str(match_contribution) + '\t' + str(region_file.split('.')[0]))



def get_allele_genotype(genotype, ref, alt):
    base_genotype = []
    # heterozygous
    if genotype == '0/1' or genotype == '1/0' or genotype == '0|1' or genotype == '1|0':
        base_genotype.append(ref)
        base_genotype.append(alt)
    # homozygous alternate
    elif genotype == '1/1' or '1|1':
        base_genotype.append(alt)
        base_genotype.append(alt)

    return base_genotype

def find_ref_matches(pileup_base_read):
    match_ref_count = 0
    for a in pileup_base_read:
        if a == '.' or a == ',':
            match_ref_count += 1
    return match_ref_count

def find_alt_matches(pileup_base_read, allele_genotype):
    match_alt_count = 0
    alt = allele_genotype[1]
    for a in pileup_base_read:
        if a.upper() == alt:
            match_alt_count += 1
    return match_alt_count


main()
#test()
