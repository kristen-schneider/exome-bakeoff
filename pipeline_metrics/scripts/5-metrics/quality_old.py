import numpy as np
import math

def main_quality(pileup_row):
    pileup_chrm = pileup_row[0]
    pileup_start = pileup_row[1]
    pileup_end = pileup_row[2]
    pileup_quality = pileup_row[6]

    return calculate_quality(pileup_quality)

def calculate_quality(pileup_quality):
    row_quality = []
    for q in pileup_quality:
        ascii_quality = float(ord(q) - ord('!'))
        phred_score = 10. ** (-ascii_quality/10.)
        quality_again = -10 * math.log10(phred_score)
        row_quality.append(ascii_quality)
    return(np.average(row_quality))
