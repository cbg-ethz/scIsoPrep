###split fltnc file per cell

import pandas as pd
import pysam
from pathlib import Path
import os
import sys
import glob
import datetime
import argparse

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq)
    letters = [complement[base] for base in bases] 
    letters = ''.join(letters)
    reverse = letters[::-1]
    return reverse

def get_bc(bc_dir, sample):
    bc_dict = {}
    for bcfile in glob.glob(bc_dir+sample+'*.txt'): #sample.whatever.txt
        df_bc = pd.read_csv(bcfile, sep = '\t')
    return list(df_bc.barcodes)

def read_fltnc(sample,bc):
    
    dic_bam_per_cell={reverse_complement(i) : [] for i in bc}

    header = []           
    Head = True
    wrong_bc = 0
    zmw = 1

    for bamfile in glob.glob(sample + "/polyA_trimming/*.bam"):
        samfile = pysam.view('-h',bamfile)
        for line in samfile.split('\n'):
            if line:
                if line[0]=='@':
                    if Head:
                        header.append(line+'\n')
                        continue
                    else: continue

                READ=line.split()
                READ[0] = '/'.join([READ[0].split('/')[0], str(zmw), READ[0].split('/')[2]])
                BC = READ[13].split(':')[2]
                if str(BC) in dic_bam_per_cell:
                    dic_bam_per_cell[BC].append(READ)
                else:
                    wrong_bc += 1
                zmw += 1
         
        Head = False
    print(zmw)
    print(wrong_bc)

    dic_bam_per_cell = {k: v for k, v in dic_bam_per_cell.items() if len(v)!=0}
    
    return dic_bam_per_cell, header
    
                
def main(args):
    
    bc = get_bc(args.bc_dir, args.sample)

    dic_bam_per_cell, header = read_fltnc(args.sample,bc)   
    
    for BC in dic_bam_per_cell:
        Path(args.sample +"/cells/" + BC+ '_' + args.sample).mkdir(parents=True, exist_ok=True)
        with open(args.sample+ "/cells/" + BC + '_' + args.sample +'/fltnc.sam','w') as sam:
            sam.writelines(header)
            for read in dic_bam_per_cell[BC]:
                sam.writelines(['\t'.join(read)+'\n'])
        pysam.samtools.view(
                "-b", "-o{}".format(args.sample + "/cells/" + BC + '_' + args.sample +'/fltnc.bam'), 
                args.sample + "/cells/" + BC + '_' + args.sample + '/fltnc.sam', catch_stdout=False)

    with open(args.sample + '/splitting_done.txt', 'w') as done:
        done.write(str(datetime.datetime.now()))


def parse_args():
    parser = argparse.ArgumentParser(
        prog='split_cells_bam.py', 
        usage='python3 split_cells_bam.py --bc_dir <bc_dir> --sample <sample> ',
        description='Divides bamfiles per cell prior to UMI deduplication'
    )
    parser.add_argument(
        '--bc_dir', type=str,
        help='Absolute or relative path(s) to directory containing barcodes files sample.whatever.txt, tsv format'
    )
    parser.add_argument(
        '--sample', type=str,
        help='sample name (should not contain ".")'
    )

    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)


