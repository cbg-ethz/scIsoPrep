import pandas as pd
import glob
import argparse
from collections import Counter
import matplotlib
#matplotlib.use("TKAgg")
from matplotlib import pyplot as plt
import numpy as np


def autopct_format_Million(values):
    def my_format(pct):
        total = sum(values)
        val = pct*total/100/1000000
        return '{:.1f}%({v:.1f}M)'.format(pct,v=val)
    return my_format

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{:.1f}%({v:,})'.format(pct,v=val)
    return my_format

def get_allInfo(allInfo):
    colnames = ['Mol', 'Gene', 'CellType', 'Barcode', 'UMI', 'PB', 'Category', 'Exon']
    allinfo = pd.read_csv(allInfo,names=colnames, header=None, sep = '\t')
    dict_count = Counter(allinfo.PB)
    filtered_PBs = [k for k in dict_count if dict_count[k]>=3]
    allinfo = allinfo[allinfo['PB'].isin(filtered_PBs)]
    return allinfo

def main(args):

    allinfo = get_allInfo(args.AllInfo)

    #piechart all molecules

    all_mol = Counter(allinfo['Category'])

    all_mol['Other'] = all_mol.pop('intergenic') + all_mol.pop('genic') + all_mol.pop('fusion') + all_mol.pop('antisense')
    all_mol['FSM'] = all_mol.pop('full-splice_match')
    all_mol['ISM'] = all_mol.pop('incomplete-splice_match')
    all_mol['NIC'] = all_mol.pop('novel_in_catalog')
    all_mol['NNC'] = all_mol.pop('novel_not_in_catalog')

    Labels = [k for k in all_mol.keys()]
    Data   = [float(v) for v in all_mol.values()]
    Colors = ['tab:blue','tab:orange','tab:red','tab:green','tab:purple']

    fig1, ax1 = plt.subplots()
    _, _, autopcts = ax1.pie(Data, 
                             labels = Labels, 
                             colors = Colors,
                             autopct=autopct_format_Million(Data),
                             pctdistance = 0.66, 
                             textprops={'fontsize': 14, 'weight':'bold'})
    plt.setp(autopcts, **{'color':'white', 'weight':'bold', 'fontsize':9.5})
    plt.savefig('All_Molecules_Category_piechart.png') 


    #piechart isoforms only

    dict_categs = dict(zip(allinfo['PB'],allinfo['Category']))

    cat_count = Counter(dict_categs.values())
    cat_count['Other'] = cat_count.pop('intergenic') + cat_count.pop('genic') + cat_count.pop('fusion') + cat_count.pop('antisense')
    cat_count['FSM'] = cat_count.pop('full-splice_match')
    cat_count['ISM'] = cat_count.pop('incomplete-splice_match')
    cat_count['NIC'] = cat_count.pop('novel_in_catalog')
    cat_count['NNC'] = cat_count.pop('novel_not_in_catalog')


    Labels = [k for k in cat_count.keys()]
    Data   = [float(v) for v in cat_count.values()]
    
    fig1, ax1 = plt.subplots()
    _, _, autopcts = ax1.pie(Data, 
                             labels = Labels, 
                             colors = Colors,
                             autopct=autopct_format(Data), 
                             pctdistance = 0.6,
                             textprops={'fontsize': 14, 'weight':'bold'})
    plt.setp(autopcts, **{'color':'white', 'weight':'bold', 'fontsize':10})
    plt.savefig('All_Isoforms_Category_piechart.png')

def parse_args():
    parser = argparse.ArgumentParser(
        prog='isoform_categ_piechart.py.py', 
        usage='python3 isoform_categ_piechart.py --AllInfo <AllInfo>',
        description='piechart of isoform categs and all molecules categs'
    )
    parser.add_argument(
        '--AllInfo', type=str, default='../AllInfo',
        help='Absolute or relative path(s) to AllInfo data in tsv format'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)

