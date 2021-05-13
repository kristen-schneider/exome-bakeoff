import argparse
import glob
import os

def get_args():

    parser = argparse.ArgumentParser()


    parser.add_argument('--in_dir',
                        dest='in_dir',
                        type=str,
                        help='path to directory storing hap.py output',
                        required=True)
    parser.add_argument('--sample',
                        dest='sample',
                        help='sample name wanted to compute over (e.g. GIAB, etc.). Default is all.',
                        required=False)
    parser.add_argument('--out_dir',
                        dest='out_dir',
                        type=str,
                        help='path to out dir (F1 score per tech combo)',
                        required=True)
    return parser.parse_args()    

def main():
    args = get_args()
    tech_combo_INDEL_SNP = make_INDEL_SNP_dict(args.in_dir, args.sample)
    write_tech_combo_file(tech_combo_INDEL_SNP, args.out_dir)

def write_tech_combo_file(tech_combo_INDEL_SNP, out_dir):
    for tech_combo in tech_combo_INDEL_SNP:
        INDEL_sub_dir = 'F1_INDEL_per_tech_combo/'
        SNP_sub_dir = 'F1_SNP_per_tech_combo/'
                
        out_INDEL = open(out_dir+INDEL_sub_dir+tech_combo+'_f1_score.txt', 'w')
        out_SNP = open(out_dir+SNP_sub_dir+tech_combo+'_f1_score.txt', 'w')

        INDEL_score = tech_combo_INDEL_SNP[tech_combo][0]
        out_INDEL.write(str(INDEL_score))
        SNP_score = tech_combo_INDEL_SNP[tech_combo][1]
        out_SNP.write(str(SNP_score))
        
def make_INDEL_SNP_dict(in_dir, sample):
    # process specified samples
    if sample == None: token = '.txt'
    else: token = sample

    # for all tech combos
    # tech_combo: [INDEL, SNP]
    tech_combo_dict = dict() 
    for f in glob.glob(in_dir):
        # line 1
        header = ''
        # line 2,3
        INDEL = []
        # line 4,5
        SNP = []
        if token in f:
            basename = os.path.basename(f)
            tech_combo_name = basename.split('_')[0]
            for line in open(f):
                if header == '': header = line
                else:
                    A = line.rstrip().split(',')
                    F1_score = float(A[13])
                    if len(INDEL) != 2: INDEL.append(F1_score)
                    else: SNP.append(F1_score)
            tech_combo_dict[tech_combo_name] = [sum(INDEL)/len(INDEL), sum(SNP)/len(SNP)]
    
        # did not find sample (e.g. GIAB)    
        else: continue
    
    return tech_combo_dict


if __name__ == '__main__':
    main()
