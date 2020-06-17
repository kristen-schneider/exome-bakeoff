import pysam
import os

artifact_file = open('/Users/krsc0813/exome-bakeoff/kristen_scripts/Artifacts/HighQual.bed', 'r')
metric_files_path = '/Users/krsc0813/exome-bakeoff/Analyses/quality/59_ACMG/downsampled/metric_files/'

def main():
    #artifact_tabix = pysam.TabixFile(artifact_file)
    header = None
    
    for line in artifact_file:
        A = line.rstrip().split('\t')
        if header == None: header = line
        else:
            print(A)
            chrm = int(A[0])
            start = float(A[1])
            end = float(A[2])

            query = [chrm, start, end]
            find_match(query)
            
def find_match(query):
    chrm = query[0]
    start = query[1]
    end = query[2]
    print(chrm)
    print(start)
    print(end)

    for mf in os.listdir(metric_files_path):
        if "_sorted.bed" in mf and ".tbi" not in mf:
            header = None
            current_mf = pysam.TabixFile(metric_files_path + mf)
        
            # debugging
            for row in current_mf.fetch(chrm, start, end, parser=pysam.asBed()):
                print(row)    
            
            #match = current_mf.fetch(chrm, start, end, parser=pysam.asBed())
            #print(match)
        else: continue


main()
