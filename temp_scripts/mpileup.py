# imports
import numpy as np
from scipy.stats import binom
import pysam


###

# tabixfile = pysam.TabixFile("example.gtf.gz")
#
# for gtf in tabixfile.fetch("chr1", 1000, 2000):
#     print (gtf.contig, gtf.start, gtf.end, gtf.gene_id)

###

# test input file and format explanation
# !!! to be changed into loop to handle all files
test_input_file = 'ex_input.txt'
##COL_1=SEQUENCE NAME
##COL_2=1-BASED COORDINATE
##COL_3=REFERENCE BASE
##COL_4=NUM READS @ POSITION
##COL_5=READ BASES
# . = match on positive strand
# , = match on negative strand
# mismatch signified by a changed in letter (lowercase = negative strand, uppercase = positive strand)
# ^ start of a read
# $ end of a read
# ! indicates some poor read quality (unsure yet how)
##COL_6=BASE QUALITIES
##COL_7=ALIGNMENT MAPPING QUALITIES

###########################
# Takes one input file and one gene: (chrom, start, end) and reports strand bias metrics
###########################
def main(input_file, chrom, start, end):
	f = open(input_file, 'r')
	o = open('test_output_file.txt', 'w')
	scratch1 = open('scratch1.txt', 'w')
	scratch2 = open('scratch2.txt', 'w')

	# to signify that we are still looking for chromosome
	# allows me to quit after having found it
	# !!! to be changed when searching across more than one gene
	not_found = True

	# to hold the binomial PMF results for every coordinate
	binomial_pmf_match = []
	binomial_pmf_mismatch = []

	# to hold the sb results for every coordinate
	sb_all = []
	gatksb_all = []

	for line in f:
		# still looking for chromosome containing gene
		if not_found:

			line = line.split('\t')

			# pull out and label useful columns
			chromosome = int(line[0])
			coordinate = int(line[1])
			referenceBase = line[2]
			numReads = int(line[3])
			reads = line[4]

			# count positive, negative-strand matches AND positive, negative-strand mismatches
			# count matches
			p_match_count = 0
			n_match_count = 0
			# count mismatches
			p_mismatch_count = 0
			n_mismatch_count = 0

			# only want to look at the specified gene
			if chromosome == chrom:
				if coordinate >= start and coordinate <= end:
					for i in range(len(reads)):

						# match + strand
						if reads[i] == '.':
							p_match_count += 1
						# match - strand
						elif reads[i] == ',':
							n_match_count += 1
						# mismatch + strand
						elif (reads[i] == 'A' or reads[i] == 'C' or reads[i] == 'T' or reads[i] == 'G'):
							p_mismatch_count += 1
						# mismatch - strand
						elif (reads[i] == 'a' or reads[i] == 'c' or reads[i] == 't' or reads[i] == 'g'):
							n_mismatch_count += 1
						# else nothing worth reporting (for now)

					# calculate binomial pmf over each coordinate
					match_bionomial = binomialPMF(p_match_count, reads, 0.5)
					mismatch_bionomial = binomialPMF(p_mismatch_count, reads, 0.5)

					# append all binomial pmf metrics to a cumulative list
					binomial_pmf_match.append(match_bionomial)
					binomial_pmf_mismatch.append(mismatch_bionomial)

					# calculate strand bias
					SB = sb(p_match_count, n_match_count, p_mismatch_count, n_mismatch_count)
					# calculate GATK strand bias
					GATKSB = gatksb(p_match_count, n_match_count, p_mismatch_count, n_mismatch_count)

					# append all sb metrics to a cumulative list
					if SB != -1: sb_all.append(SB)
					if GATKSB != -1: gatksb_all.append(GATKSB)


					# debug
					print >> scratch1, p_match_count, n_match_count, p_mismatch_count, n_mismatch_count, len(reads)

			# found chromosome containing gene, exit loop
			elif chromosome > chrom: not_found = False

	print >> scratch2, binomial_pmf_match, '\n', binomial_pmf_mismatch

	# calculate summary metrics for all the coordinates in one gene
	binomial_pmf_match_metrics = getSummaryMetrics(binomial_pmf_match)
	binomial_pmf_mismatch_metrics = getSummaryMetrics(binomial_pmf_mismatch)
	sb_metrics = getSummaryMetrics(sb_all)
	gatksb_metrics = getSummaryMetrics(gatksb_all)

	# binomial pmf
	print >> o, "\nMATCH", '\n', "max: ", binomial_pmf_match_metrics[0], '\n', "min: ", binomial_pmf_match_metrics[1], '\n', "mean: ", \
		binomial_pmf_match_metrics[2], '\n', "median: ", binomial_pmf_match_metrics[3], '\n', "std: ", binomial_pmf_match_metrics[4]
	print >> o, "\n\nMISMATCH", '\n', "max: ", binomial_pmf_mismatch_metrics[0], '\n', "min: ", binomial_pmf_mismatch_metrics[1], '\n', "mean: ",\
		binomial_pmf_mismatch_metrics[2], '\n', "median: ", binomial_pmf_mismatch_metrics[3], '\n', "std: ", binomial_pmf_mismatch_metrics[4]
	# sb
	print >> o, "\n\nSB", '\n', "max: ", sb_metrics[0], '\n', "min: ", sb_metrics[1], '\n', "mean: ", sb_metrics[2],\
		'\n', "median: ", sb_metrics[3], '\n', "std: ", sb_metrics[4]
	print >> o, "\n\nGATK-SB", '\n', "max: ", gatksb_metrics[0], '\n', "min: ", gatksb_metrics[1], '\n', "mean: ", gatksb_metrics[2], \
		'\n', "median: ", gatksb_metrics[3], '\n', "std: ", gatksb_metrics[4]

	f.close()

###########################
# Calculate strand bias given forward and reverse match and mismatch counts
# this is a modification of strand bias as discussed in:
# 'The effect of strand bias in Illumina short-read sequencing data' (Y. Guo, et. al.)
# The modification comes in because I do not define a "minor" allele so I am counting all mismatches
###########################
def sb(p_match, n_match, p_mismatch, n_mismatch):
	a = float(p_match)
	b = float(p_mismatch)
	c = float(n_match)
	d = float(n_mismatch)

	try: SB = abs((b/(a+b))-(d/(c+d)))/((b+d)/(a+b+c+d))
	except ZeroDivisionError: return -1

	return SB


###########################
# Calculate GATK strand bias given forward and reverse match and mismatch counts
# this is a modification of GATK strand bias as discussed in:
# 'The effect of strand bias in Illumina short-read sequencing data' (Y. Guo, et. al.)
# The modification comes in because I do not define a "minor" allele so I am counting all mismatches
###########################
def gatksb(p_match, n_match, p_mismatch, n_mismatch):
	a = float(p_match)
	b = float(p_mismatch)
	c = float(n_match)
	d = float(n_mismatch)
	try: GATKSB = max(((b/(a+b))*(c/(c+d)))/((a+c)/(a+b+c+d)), ((d/(c+d))*(a/(a+b)))/((a+c)/(a+b+c+d)))
	except ZeroDivisionError: return -1

	return GATKSB


###########################
# Calculates the probability mass function of the binomial distribution given
# positive strand counts (successes), Read Bases (col 5), and a probability of success (0.5)
###########################
def binomialPMF(p_count, reads, prob_success):

	# binomial distribution: what is the probability of x "successes" given n trials with a success probability of p?
	# p_count: x = number of successes = counts on positive strand
	# num_reads: n = number of trials = number of reads
	# prob_success: p = probability of success = 0.5

	x = p_count
	n = len(reads)
	p = prob_success

	return binom.pmf(x, n, p)

###########################
# Reports summary metrics of the binomial pmf
###########################
def getSummaryMetrics(array_count):

	# max, min, mean, median, mode
	metrics = []

	max = np.max(array_count)
	min = np.min(array_count)
	mean = np.average(array_count)
	median = np.median(array_count)
	std = np.std(array_count)
	# mode_counts = np.bincount(array_count)
	# mode = np.argmax(mode_counts)

	metrics.extend([max, min, mean, median, std])

	return metrics

###########################
# Retrieves the fifth column using pysam
###########################

tbx = pysam.TabixFile('AgilentQXT-IDT-0720ME25_S12_L001.bqsr.bam.bed.gz')

for row in tbx.fetch('1', 10000, 10100, parser=pysam.asBed()):
	print row[5]

# testing with REEP1 Duplication given by Rebecca
main(test_input_file, 2, 86305963, 86510812)

## garbage testing
# kristen = [2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 3.0, 3.0, 3.0]
# print getSummaryMetrics(kristen)
#
# print binomialPMF(3, [0,0,0,0,0,0,7], 0.5)
print sb(10, 9, 3, 4)
print gatksb(10, 9, 3, 4)


