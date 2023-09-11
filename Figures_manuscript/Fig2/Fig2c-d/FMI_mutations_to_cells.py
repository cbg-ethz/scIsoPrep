from collections import defaultdict
import pandas as pd
import argparse
import pysam
import re


def get_cigartuple(CIGAR):
    l=re.split('(\d+)', CIGAR)
    cigartuples=[[int(l[i+1]),l[i+2]] for i in range(0,len(l)-2,2)]
    return cigartuples


def get_readpos(cigartuples,pos_read,pos_genome,mutpos):
    
    for i in cigartuples:
        if i[1]=='H':
            continue
        elif i[1]=='S':
            pos_read += i[0]
        elif i[1]=='M':
            pos_genome += i[0] #progress on genome
            if pos_genome >= mut_pos: #if exon exceed mutation pos, find mutation pos on read
                pos_read += i[0] - (pos_genome - mut_pos) #mutation pos on read = actual pos + distance to pos in exon 
                return pos_read
            pos_read += i[0]    
        elif i[1]=='N':
            pos_genome += i[0]
            if pos_genome >= mut_pos:
                return -1
        elif i[1]=='I':
            pos_read += i[0]
        elif i[1]=='D':
            pos_genome += i[0]
            if pos_genome >= mut_pos:
                return -1

def get_ReadBase(READ,mut_pos):

    SEQ = READ.query_sequence
    START = READ.reference_start
    CIGAR = READ.cigarstring
    NAME = READ.query_name
    
    CELL = '_'.join(NAME.split('_')[1:4])

    pos_genome = START #initialize position on genome
    pos_read = 0 #initialize position on read
    
    cigartuples = get_cigartuple(CIGAR)
    
    pos_read = get_readpos(cigartuples,pos_read,pos_genome,mutpos)
    
    if pos_read == -1:
        return 'X'
    
    base = SEQ[pos_read]
    return CELL, base


def main(args):
    
    muts = pd.read_csv(args.csv)
    sam = pysam.AlignmentFile(args.bam, "rb")
    MutatedCells = defaultdict(lambda: [])
    
    for mut in muts.iterrows():
        mut_base = mut['Mut']
        mchr, mut_pos = mut['Position'].split(':')
        
        for READ in sam.fetch(mchr,mut_pos,mut_pos+1):
            
            CELL, base = get_ReadBase(READ,mut_pos)
            if base == mut_base:
                MutatedCells[CELL].append(base)
                
                
    with open(args.out) as f:
        f.write('Cell,Mutations\n')
        for cell in MutatedCells:
            mutationChain = ';'.join(MutatedCells[cell])
            f.write('{},{}'.format(cell,mutationChain))
                
def parse_args():
    parser = argparse.ArgumentParser(
        prog='FMI_mutations_to_cells.py', 
        usage='python3 FMI_mutations_to_cells.py --bam <sorted.bam> --csv FMI_mutations.csv',
        description='Creates csv files of mutated cell names and their mutations'
    )
    parser.add_argument(
        '--csv', type=str, default='FMI_mutations.csv',
        help='Absolute or relative path(s) to AllInfo data in tsv format'
    )
    parser.add_argument(
        '--bam', type=str, default='merged_all.sorted_by_name.bam',
        help='Absolute or relative path(s) to bam file containing all reads'
    )
    parser.add_argument(
        '-o','--out', type=str, default='MutatedCells.csv',
        help='output csv file listing cells mutated'
    )


    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
