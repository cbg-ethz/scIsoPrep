### merge sqanti_filtered of all cells with barcode
import pandas as pd
import os
import argparse
import glob 


def get_filtered_PBs(cell_folder):
    f = os.path.join(cell_folder, "sqanti3_qc/x.collapsed_classification.filtered_lite_classification.txt")
    df_PB = pd.read_csv(f, sep = '\t')
    PBs = list(df_PB.isoform)
    return set(PBs)
    
def group_to_mol(cell_folder):
    f = os.path.join(cell_folder, "mapping/x.collapsed.group.txt")
    df_grp = pd.read_csv(f, sep = '\t', header = None)

    PB_to_mol = dict(zip(df_grp.iloc[:, 0], df_grp.iloc[:, 1])) #zip dict PB to associated molecules
    PB_to_mol = {key: list(map(str, value.split(','))) for key, value in PB_to_mol.items()}
    return PB_to_mol


def PB_to_filtered_mol(PBs, PB_to_mol):
    final_mol = []
    
    for PB in PB_to_mol:
        if PB in PBs:
            for mol in PB_to_mol[PB]:
                final_mol.append(mol)
    
    return set(final_mol)

def dedup_info_csv(cell_folder, sample, dedup_info_list, bc):

    f = os.path.join(cell_folder,"dedup/dedup.info.csv")
    dedup_info = pd.read_csv(f, sep = '\t')
    dedup_info['id'] = '_'.join([str(dedup_info['id']), bc, sample])

    dedup_info_list.append(dedup_info)
    
    return dedup_info_list
        

def collect_sam(cell_folder, sample, filtered_mol, bc):

    collected_sam = []
    with open(cell_folder + "mapping/mapped.fasta.sam", 'r') as sam:
        reads = sam.readlines()

        for read in reads:
            s = read.split('\t')
            if s[0] in filtered_mol:
                s[0] = '_'.join([s[0],bc,sample])
                collected_sam.append('\t'.join(s))
            
    return collected_sam


def get_header(cell_folder):
    header = []
    f = os.path.join(cell_folder,"mapping/mapped.fasta.sam")
    with open(f, 'r') as sam:
        reads = sam.readlines()

        for read in reads:
            if read[0]=='@':
                header.append(read)
            else:
                break
    return header
            
def main(args):
    sample = str(args.sample)
    final_sam = []
    dedup_info_list = []
    
    for cell_folder in glob.glob("{}/cells/*/".format(sample)):
        bc = cell_folder.split('/')[-2]
        dedup_info_list = dedup_info_csv(cell_folder, sample, dedup_info_list, bc)
        PB = get_filtered_PBs(cell_folder)
        PB_to_mol = group_to_mol(cell_folder)
        filtered_mol = PB_to_filtered_mol(PB, PB_to_mol)
        collected_final_sam = collect_sam(cell_folder, sample, filtered_mol, bc)
        final_sam += collected_final_sam
    header = get_header(cell_folder)
    
    final_dedup = pd.concat(dedup_info_list)

    final_dedup.to_csv("{}/mapping/{}.dedup.info.csv".format(sample,sample), sep = '\t', index=False)
                             
    with open("{}/mapping/{}.merge.sam".format(sample,sample), 'w') as final:
        for read in header:                      
            final.write(read)                    
        for read in final_sam:
            try:
                final.write(read)
            except TypeError:
                print(read)
                
                                 
    return

def parse_args():
    parser = argparse.ArgumentParser(
        prog='merge_cells_bc.py', 
        usage='python3 merge_cells_bc.py --sample <sample>',
        description='Merges files from all cells of a sample'
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




