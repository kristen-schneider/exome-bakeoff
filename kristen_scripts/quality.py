import numpy as np

def main_quality(pileup_row):
    pileup_chrm = pileup_row[0]
    pileup_start = pileup_row[1]
    pileup_end = pileup_row[2]
    pileup_quality = pileup_row[6]

    return calculate_quality(pileup_quality)

def calculate_quality(pileup_quality):
    row_quality = []
    for q in pileup_quality:
        row_quality.append(ord(q) - ord('!'))
    return(np.average(row_quality))
