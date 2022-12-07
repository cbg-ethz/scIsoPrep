import pandas as pd
import sys
import argparse


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq)
    letters = [complement[base] for base in bases] 
    letters = ''.join(letters)
    reverse = letters[::-1]
    return reverse

def get_cell_type(ctypes):
    all_ctypes = {}
    all_bc_rev = {}
    for c in ctypes:
        sample = c.split('/')[-1].split('.')[0]
        df = pd.read_csv(c, sep = '\t')
        bc_rev = {reverse_complement(bc) : bc for bc in df['barcodes']}
        cell_types = dict(zip(df['barcodes'], sample + '_' + df['celltype_final']))
        all_ctypes.update(cell_types)
        all_bc_rev.update(bc_rev)
    return all_ctypes, all_bc_rev

def main(args):

    dico_ctypes,dico_bc = get_cell_type(args.cell_types)

    dico_class = pd.read_csv(args.categ, sep='\t').set_index('isoform').to_dict()
    
    df_group = pd.read_csv(args.grp, names=['pb','mol'], sep='\t')
    
    AllInfo = {}
    
    for index, row in df_group.iterrows():
        pb = row['pb']
        try:
            a = dico_class['associated_gene'][pb]
        except KeyError:
            continue
        for mol in row['mol'].split(','):
            rev_bc = mol.split('_')[1]
            bc = dico_bc[rev_bc]
            umi, sep, ID = mol.partition('_')
            
            AllInfo[mol] = {'mol': mol, 
                            'gene': dico_class['associated_gene'][pb],
                            'sample': dico_ctypes[bc],
                            'bc': bc,
                            'umi': umi,
                            'isoform': pb,
                            'category': dico_class['structural_category'][pb],
                            'exon': dico_class['exons'][pb] 
                           }
    
    df_info = pd.DataFrame.from_dict(AllInfo, orient='index')
    
    df_info.to_csv(args.out,
        sep = '\t', index=False, header=False)
    
    
def parse_args():
    parser = argparse.ArgumentParser(
        prog='AllInfo_V2.py', 
        usage='python3 AllInfo_V2.py --grp <sqanti_group.gff> --csv <coldata.csv> -o <OUT> [options]',
        description='Creates AllInfo file from SQANTI output for scisorseq differential isoform expression'
    )
    parser.add_argument(
        '--grp', type=str,
        help='Absolute or relative path(s) to SQANTI x.collapsed.group.txt data'
    )
    parser.add_argument(
        '--cell_types', type=str, nargs='+',
        help='Absolute or relative path(s) to coldata.csv data in CSV format'
    )
    parser.add_argument(
        '-o', '--out', type=str, default='AllInfo',
        help='Output file name. Default = AllInfo'
    )
    parser.add_argument(
        '--categ', type=str,
        help='Absolute or relative path(s) to SQANTI x.collapsed_classification.filtered_lite_classification.txt data'
    )
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
