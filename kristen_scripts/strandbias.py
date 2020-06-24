def main_strandbias(pileup_row):
    pileup_chrm = pileup_row[0]
    pileup_start = pileup_row[1]
    pileup_end = pileup_row[2]
    pileup_sb = pileup_row[5]
    counts = get_counts(pileup_sb)
    
    return calculate_strandbias(counts)

def get_counts(pileup_sb):
    f_match = 0
    f_mismatch = 0
    r_match = 0
    r_mismatch = 0
    for i in pileup_sb:
        if i == '.': f_match+=1
        elif i == ',': r_match+=1
        elif i in 'ACTGN': f_mismatch+=1
        elif i in 'actgn': r_mismatch+=1
        else: continue
    return [f_match, f_mismatch, r_match, r_mismatch]


def calculate_strandbias(counts):
    a = float(counts[0])
    b = float(counts[1])
    c = float(counts[2])
    d = float(counts[3])
    
    try: sb = (c+d)/(a+b+c+d)
    except ZeroDivisionError:
        return -1
    
    #strand_bias = []
    #for c_i in counts:
    #    curr_position = c_i[0]
    #    curr_counts = c_i[1]
    #    #print (c_i)
    #    a = float(curr_counts[0])
    #    b = float(curr_counts[1])
    #    c = float(curr_counts[2])
    #    d = float(curr_counts[3])

    #    try: sb = (c+d)/(a+b+c+d)
    #    except ZeroDivisionError:
    #        return -1
    #    curr_row = [curr_position, sb]
    #    strand_bias.append(curr_row)

    return sb
