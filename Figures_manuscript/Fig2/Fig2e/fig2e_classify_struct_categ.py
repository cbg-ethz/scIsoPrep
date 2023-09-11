from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import pandas as pd
import argparse

Colors = ['tab:orange','tab:green','tab:purple']

def get_allInfo(allInfo):
    colnames = ['Mol', 'Gene', 'CellType', 'Barcode', 'UMI', 'PB', 'Category', 'Exon']
    allinfo = pd.read_csv(allInfo,names=colnames, header=None, sep = '\t')
    allinfo['Barcode'] = allinfo['Barcode'] + '_' + allinfo['CellType'].str.split('_',expand=True,)[0]  + '_' + allinfo['CellType'].str.split('_',expand=True,)[1]
    allinfo['CellType'] = allinfo['CellType'].str.split('_',expand=True,)[2]
    allinfo = allinfo[allinfo['Category'].isin(
        ['novel_not_in_catalog', 'novel_in_catalog', 'full-splice_match'])]
    return allinfo

def stacked_barplot_normalized(df,df_labels,ax1,xlab):

    #getting raw value labels
    a=[]
    for c in df_labels.loc[:, df_labels.columns != 'categ']:
        a+=list(df_labels[c])
        
    ax = df.plot(
        x = 'categ',
        kind = 'barh',
        stacked = True,
        color = Colors,
        mark_right = True,
        fontsize=20,
    ax = ax1)

    i=0
    for rect in ax.patches:
        # Find where everything is located
        height = rect.get_height()
        width = rect.get_width()
        x = rect.get_x()
        y = rect.get_y()

        # The height of the bar is the data value and can be used as the label
        #label_text = f'{width}'  # f'{height:.2f}' to format decimal values
        label_text = a[i]
        i+=1
        # ax.text(x, y, text)
        label_x = x + width / 2
        label_y = y + height / 2

        # plot only when height is greater than specified value
        if height > 0:
            ax.text(label_x, label_y, label_text, color='white', weight='bold', ha='center', va='center', fontsize=15)

    ax.legend(bbox_to_anchor=(1.05, 0.1), loc='upper left', borderaxespad=0.,fontsize=20)
    ax.get_legend().remove()

    ax.set_ylabel('', fontsize=0)
    #ax.tick_params(axis='y', which='both',left=False, right=False, labelleft=False)
    ax.tick_params(axis='x', which='major', labelsize=16)
    ax.set_xlabel(xlab, fontsize=16)
    ax.invert_yaxis()


def set_up_df(count_categ):

    df = pd.DataFrame(count_categ)
    df = df.rename(columns={'full-splice_match':'FSM', 'novel_in_catalog':'NIC', 'novel_not_in_catalog': 'NNC'})
    df = df.loc[['HGSOC','Mesothelial.cells','Myeloid.cells','T.NK.cells','Fibroblasts', 'B.cells','Endothelial.cells']]

    sum_tot = df.sum(axis=1)

    dfb = df.div(sum_tot, axis=0).multiply(100)
    df.index.name = 'categ'
    df.reset_index(inplace=True)
    dfb.index.name = 'categ'
    dfb.reset_index(inplace=True)

    df = df.reindex(sorted(df.columns), axis=1)
    dfb = dfb.reindex(sorted(dfb.columns), axis=1)

    return df, dfb, sum_tot


def main(args):

    allinfo = get_allInfo(args.AllInfo)

    iso_ctype = defaultdict(lambda: {'ctypes':set(),'mols':0})
    count_categ = defaultdict(lambda: defaultdict(lambda:set()))

    for row in allinfo.itertuples(index = False):
        pb = row[5]
        ctype = row[2]
        categ = row[6]
        if ctype == 'uncertain':
                continue
        iso_ctype[pb]['ctypes'].add(ctype)
        iso_ctype[pb]['mols']+=1
        iso_ctype[pb]['categ'] = categ

    for pb in iso_ctype:
        nmols = iso_ctype[pb]['mols']
        if nmols < args.filter:
            continue
        categ = iso_ctype[pb]['categ']
        for ctype in iso_ctype[pb]['ctypes']:
            count_categ[categ][ctype].add(pb)
            

    for categ in count_categ:
        for ctype in count_categ[categ]:
            count_categ[categ][ctype] = len(count_categ[categ][ctype])

    df, dfb, sum_tot = set_up_df(count_categ)

    fig, axes = plt.subplots(figsize=(7, 5), dpi = 300)
    stacked_barplot_normalized(dfb,df,axes,"Percentage of structural categories")
    fig.tight_layout()   
    fig.savefig(args.out,bbox_inches = 'tight')


def parse_args():
    parser = argparse.ArgumentParser(
        prog='fig2e_classify_struct_categ.py', 
        usage='python3 classifier.py --AllInfo <AllInfo> -o <out.png> [options]',
        description='Takes AllInfo and SQANTI classification to classify celltype specific isoforms'
    )
    parser.add_argument(
        '--AllInfo', type=str, default = '../AllInfo',
        help='Absolute or relative path(s) to AllInfo data in tsv format'
    )
    parser.add_argument(
        '-o', '--out', type=str, default='Fig3e_isoforms_overall.png',
        help='Output directory'
    )

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
