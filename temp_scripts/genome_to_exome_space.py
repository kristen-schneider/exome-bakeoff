import sys

def geno_to_exo_main(regions_file, old_plot_file, output):
    print('hello')
    in_file = open(old_plot_file, 'r')
    out_file = open(output, 'a')
    out_file.truncate(0)

    regions = []
    for l in open(regions_file):
        A = l.rstrip().split()
        chrom = A[0]
        start = int(A[1])
        end = int(A[2])
        regions.append([chrom,start,end])

    tally = 0
    for r in regions:
        r.append(tally)
        tally += r[2] - r[1]
    
    for l in in_file:
        A = l.rstrip().split()
        start = int(A[1])
        end = int(A[2])
        score = float(A[3])
        i = bsearch(start, regions, -1, len(regions))
        offset = regions[i][3] + (start - regions[i][1])
        print (offset, score, file=out_file)



def bsearch(key, D, lo, hi):

    i = 0
    while (hi - lo > 1):
        i+=1
        mid = (hi + lo) // 2
        #print mid, key, D[mid], key < D[mid][1]
        # D[mid] -> [chrm, start, end, offset]
        if ( key >= D[mid][1] and key <= D[mid][2] ):
            return mid
        if ( key < D[mid][1] ):
            hi = mid
        else:
            assert( key > D[mid][2] )
            lo = mid
    return hi

#regions_file = sys.argv[1]
#pile_up_file = sys.argv[2]

