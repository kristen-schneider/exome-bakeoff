import os

artifact_csv = open('/Users/krsc0813/exome-bakeoff/kristen_scripts/Artifacts/HighQual.csv', 'r')

def main():
    artifact_bed = open('/Users/krsc0813/exome-bakeoff/kristen_scripts/Artifacts/HighQual.bed', 'a')
    artifact_bed.truncate(0)
    header = None
    
    for line in artifact_csv:
        if header == None: header = line
        else:
            A = line.rstrip().split(',')
            artifact_bed.write(A[0] + '\t' + A[1] + '\t' + str(int(A[1])+1) + '\n')
    
        

main()
