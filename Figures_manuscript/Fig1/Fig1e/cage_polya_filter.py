import pandas as pd
import glob
import argparse
from collections import Counter
import matplotlib
#matplotlib.use("TKAgg")
from matplotlib import pyplot as plt
import numpy as np

#generates polyA/cage support barchart (Fig1e)
#GENCODE FL / MANE gff files were generated by filtering isoforms in gencode_v36
#Input bed_polyA files were generated using:
# bedtools closest -a data.gff -b atlas.clusters.2.0.GRCh38.96.bed
#Input bed_polyA files were generated using:
# bedtools closest -a data.gff -b atlas.clusters.2.0.GRCh38.96.bed
#GENCODE FL / MANE were generated by filtering isoforms classified as such in gencode_v36

def get_allInfo(allInfo):
    colnames = ['Mol', 'Gene', 'CellType', 'Barcode', 'UMI', 'PB', 'Category', 'Exon']
    allinfo = pd.read_csv(allInfo,names=colnames, header=None, sep = '\t')
    dict_count = Counter(allinfo.PB)
    filtered_PBs = [k for k in dict_count if dict_count[k]>=3]
    allinfo = allinfo[allinfo['PB'].isin(filtered_PBs)]
    dict_categ = dict(zip(allinfo.PB, allinfo.Category))
    return allinfo, dict_categ, filtered_PBs

def get_bed_gencode_polya(bed_file):
    bed = pd.read_csv(bed_file, names = list(range(21)), header=None, sep = '\t')
    iso_types = ['protein_coding', 'nonsense_mediated_decay', 'retained_intron', 'processed_transcript', 'processed_pseudogene','transcribed_processed_pseudogene','transcribed_unprocessed_pseudogene','transcribed_unitary_pseudogene']
    bed.columns = bed.columns.astype(str)
    bed['iso_type'] = [i.split(';')[4].split('"')[1] for i in list(bed['8'])]
    bed = bed[bed['iso_type'].isin(iso_types)]
    PBs = [i.split(';')[1].split('"')[1] for i in list(bed['8'])]
    bed['8'] = PBs
    #bed_filter = bed.sort_values('18').drop_duplicates('8').sort_index()
    #bed_filter = bed_filter[bed_filter['18']<=50]
    dict_dist = dict(zip(bed['8'],bed['20']))
    tot = list(dict_dist.keys())
    qualif = [k for k in dict_dist if dict_dist[k]==0]

    return qualif, tot

def get_bed_gencode_cage(bed_file, dist):
    bed = pd.read_csv(bed_file, names = list(range(19)), header=None, sep = '\t')
    bed.columns = bed.columns.astype(str)                                                                                                                
    PBs = [i.split(';')[1].split('"')[1] for i in list(bed['8'])]
    bed['8'] = PBs
    rows = []
    for row in bed.itertuples(index=False):
        if row[6] == '+':
            a = row[3]
        elif row[6] == '-':
            a = row[4]
        row = list(row)
        if abs(row[10] - a) <=dist or abs(row[11] - a) <=dist:
            #row[18] = min(abs(row[10] - a),abs(row[11] - a))
            row[18] = 0
        else:
            row[18] = dist+1
        rows.append(row)
    bed_final = pd.DataFrame(rows, columns =list(range(19)))
    bed_final.columns = bed_final.columns.astype(str)
    bed_final = bed_final.sort_values('18').drop_duplicates('8').sort_index()

    dict_dist = dict(zip(bed_final['8'],bed_final['18']))
    tot = list(dict_dist.keys())

    qualif = [k for k in dict_dist if dict_dist[k]<=dist]

    return qualif, tot

def get_bed_cage(bed_file, dist):
    bed = pd.read_csv(bed_file, names = list(range(19)), header=None, sep = '\t')
    bed.columns = bed.columns.astype(str)                                                                                                                
    PBs = [i.split(';')[0].split('"')[1] for i in list(bed['8'])]
    bed['8'] = PBs
    #bed_filter = bed.sort_values('18').drop_duplicates('8').sort_index()
    #bed_filter = bed_filter[bed_filter['18']<=50]
    rows = []
    for row in bed.itertuples(index=False):
        if row[6] == '+':
            a = row[3]
        elif row[6] == '-':
            a = row[4]
        row = list(row)
        if abs(row[10] - a) <=dist or abs(row[11] - a) <=dist:
            #row[18] = min(abs(row[10] - a),abs(row[11] - a))
            row[18] = 0
        else:
            row[18] = dist+1
        rows.append(row)
    bed_final = pd.DataFrame(rows, columns =list(range(19)))
    bed_final.columns = bed_final.columns.astype(str)
    bed_final = bed_final.sort_values('18').drop_duplicates('8').sort_index()

    dict_dist = dict(zip(bed_final['8'],bed_final['18']))
    
    return dict_dist

def get_bed_polyA(bed_file):
    bed = pd.read_csv(bed_file, names = list(range(21)), header=None, sep = '\t')
    bed.columns = bed.columns.astype(str)                                                                                                                
    PBs = [i.split(';')[0].split('"')[1] for i in list(bed['8'])]
    bed['8'] = PBs
    #bed_filter = bed.sort_values('18').drop_duplicates('8').sort_index()
    #bed_filter = bed_filter[bed_filter['18']<=50]
    dict_dist = dict(zip(bed['8'],bed['20']))
    return dict_dist

def get_classif(classif):
    df_classif = pd.read_csv(classif, sep = '\t')
    TSS = df_classif[abs(df_classif['diff_to_gene_TSS']) <=100]['isoform']
    df_classif = df_classif[df_classif['perc_A_downstream_TTS'] < 80.0]
    df_classif = df_classif[abs(df_classif['diff_to_gene_TTS']) <=100]
    TTS = df_classif['isoform']
    return TSS,TTS

def compare(qual1,tot1,qual2,tot2):
   qual =  [x for x in qual1 if x in qual2]
   tot = list(set(tot1+tot2))
   return len(qual)/len(tot)*100


def main(args):
    
    percentages = {}

    TSS_PB, TTS_PB = get_classif(args.classif) 

    dist_cage = get_bed_cage(args.PBbed_cage, arg.dist_cage)
    cage_PB = [PB for PB in dist_cage if dist_cage[PB] <= args.dist] 
    keep_CAGE = set(TSS_PB + cage_PB)

    dist_polyA = get_bed_polya(args.PBbed_polya)
    polyA_PB = [PB for PB in dist_polyA]
    keep_polyA = set(TTS_PB + polyA_PB)

    keep = keep_CAGE
    keep.update(keep_polyA)


    qualif_gencode_polya, tot_gencode_polya = get_bed_gencode_polya(args.genbed_polyA)
    qualif_gencode_cage, tot_gencode_cage = get_bed_gencode_cage(args.genbed_cage, args.dist)
    percentages['gencode.full'] = [compare(qualif_gencode_polya, tot_gencode_polya, qualif_gencode_cage, tot_gencode_cage)]

    qualif_gencode_FL_polya, tot_gencode_FL_polya = get_bed_gencode_polya(args.FLbed_polyA)
    qualif_gencode_FL_cage, tot_gencode_FL_cage = get_bed_gencode_cage(args.FLbed_cage, args.dist)
    percentages['gencode.FL'] = [compare(qualif_gencode_FL_polya, tot_gencode_FL_polya, qualif_gencode_FL_cage, tot_gencode_FL_cage)]

    qualif_gencode_MANE_polya, tot_gencode_MANE_polya = get_bed_gencode_polya(args.MANEbed_polyA)
    qualif_gencode_MANE_cage, tot_gencode_MANE_cage = get_bed_gencode_cage(args.MANEbed_cage, args.dist)
    percentages['gencode.MANE'] = [compare(qualif_gencode_MANE_polya, tot_gencode_MANE_polya, qualif_gencode_MANE_cage, tot_gencode_MANE_cage)]

    colnames = ['Mol', 'Gene', 'CellType', 'Barcode', 'UMI', 'PB', 'Category', 'Exon']
    allinfo = pd.read_csv(args.AllInfo,names=colnames, header=None, sep = '\t')
    allinfo = allinfo.drop_duplicates(subset=['PB'])
    dict_categ = Counter(allinfo.Category)

    final_allinfo = allinfo[['PB'].isin(keep)]
    dict_categ_final = Counter(final_allinfo.Category)

    percentages['FSM'] = [(dict_categ_final['full-splice_match'] / dict_categ['full-splice_match'])*100]
    percentages['ISM'] = [(dict_categ_final['incomplete-splice_match'] / dict_categ['incomplete-splice_match'])*100]
    percentages['NIC'] = [(dict_categ_final['novel_in_catalog'] / dict_categ['novel_in_catalog'])*100]
    percentages['NNC'] = [(dict_categ_final['novel_not_in_catalog'] / dict_categ['novel_not_in_catalog'])*100]

    dict_categ_final['Other'] = dict_categ_final['intergenic'] + dict_categ_final['genic'] + dict_categ_final['fusion'] + dict_categ_final['antisense']
    dict_categ['Other'] = dict_categ['intergenic'] + dict_categ['genic'] + dict_categ['fusion'] + dict_categ['antisense']

    percentages['Other'] = [(dict_categ_final['Other'] / dict_categ['Other'])*100]

    df = pd.DataFrame(percentages)

    cat_plot = {}

    for column in df:
        cat_plot[column] = df[column][0]


    l = ['GENCODE.MANE', 'NIC','FSM','NNC','Other','ISM','GENCODE.FL','GENCODE.all']
    Data = []
    Labels = []
    for i in l:
        Labels.append(i)
        Data.append(float(cat_plot[i]))  

    #cat_plot = {key: value for key, value in cat_plot.items()}         
    #Labels = [k for k in cat_plot.keys()]
    Colors = ['grey','tab:green','tab:orange','tab:purple','tab:blue','tab:red','grey','grey']
    #Data   = [float(v) for v in cat_plot.values()]

    plt.figure(figsize=(4, 3), dpi = 600) 

    plt.bar(Labels, Data, color = Colors)
    plt.ylim([0,100])
    plt.ylabel("Isoforms with 5' and 3'\n support [%]", fontsize=13)
    plt.xticks(rotation=45, fontsize=11, ha = 'right', rotation_mode="anchor")
    plt.yticks(fontsize=11)
    plt.tight_layout()
    plt.savefig(args.out, bbox_inches='tight')

def parse_args():
    parser = argparse.ArgumentParser(
        prog='cage_polya_filter.py', 
        usage='python3 cage_polya_filter.py  -o <out.png> [options]',
        description='compute polyA and cage filtered '
    )
    parser.add_argument(
        '--AllInfo', type=str, default='../AllInfo',
        help='Absolute or relative path(s) to AllInfo data in tsv format'
    )
    parser.add_argument(
        '--genbed_polyA', type=str, default = './closest.polyA.gencodev36.bed',
        help='Absolute or relative path(s) to closest to gencode last exon polyA data in bed format'
    )
    parser.add_argument(
        '--MANEbed_polyA', type=str, default = './closest.polyA.MANE.gencodev36.bed',
        help='Absolute or relative path(s) to closest to gencode last exon polyA data in bed format'
    )
    parser.add_argument(
        '--FLbed_polyA', type=str, default = './closest.polyA.FL.gencodev36.bed',
        help='Absolute or relative path(s) to closest to gencode full-length last exon polyA data in bed format'
    )
    parser.add_argument(
        '--PBbed_cage', type=str, default = './closest.cage_first_exon.bed',
        help='Absolute or relative path(s) to closest to gencode first exon CAGE data in bed format'
    )
    parser.add_argument(
        '--PBbed_polyA', type=str, default = './closest.polyA_last_exon.bed',
        help='Absolute or relative path(s) to closest to gencode first exon CAGE data in bed format'
    )
    parser.add_argument(
        '--genbed_cage', type=str, default = './closest.CAGE_gencode.bed',
        help='Absolute or relative path(s) to closest to gencode first exon CAGE data in bed format'
    )
    parser.add_argument(
        '--MANEbed_cage', type=str, default = './closest.CAGE.MANE.gencode.bed',
        help='Absolute or relative path(s) to closest to gencode last exon CAGE data in bed format'
    )
    parser.add_argument(
        '--FLbed_cage', type=str, default = './closest.CAGE.FL.gencode.bed',
        help='Absolute or relative path(s) to closest to gencode full-length last exon polyA data in bed format'
    )
    parser.add_argument(
        '--dist_cage', type=int, default=100,
        help='max distance from CAGE peak'
    )
    parser.add_argument(
        '-o', '--out', type=str, default='./fig1e.png',
        help='Output filtered AllInfo file'
    )


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)




