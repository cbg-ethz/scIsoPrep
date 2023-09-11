import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
from collections import Counter


colors = ['#A034F0','#FFA501','#FF2500','#A52A2A','#0433FF'] 
labels=['HGSOC','Fibroblast','Mesothelial.cells','Myeloid.cells', 'T.cells']

def get_allInfo(allInfo):
    colnames = ['Mol', 'Gene', 'CellType', 'Barcode', 'UMI', 'PB', 'Category', 'Exon']
    allinfo = pd.read_csv(allInfo,names=colnames, header=None, sep = '\t')
    allinfo['IGF2BP2-TESPA1'] = 1
    allinfo['BarcodeTotal'] = allinfo['Barcode'] + '_' + allinfo['CellType'].str.split('_',expand=True,)[0]  + '_' + allinfo['CellType'].str.split('_',expand=True,)[1]
    allinfo['Patient'] = allinfo['CellType'].str.split('_',expand=True,)[1]
    allinfo['Patient'] = allinfo['Patient'].str.replace("B486", "Patient 1")
    allinfo['Patient'] = allinfo['Patient'].str.replace("B497", "Patient 2")
    allinfo['Patient'] = allinfo['Patient'].str.replace("B500", "Patient 3")
    allinfo['CellType'] = allinfo['CellType'].str.split('_',expand=True,)[0]
    allinfo['CellType'] = allinfo['CellType'].str.replace(".Om", "")
    return allinfo

def main(args)
    allinfo = get_allInfo(args.AllInfo)
    allinfo = allinfo[allinfo['CellType']!='uncertain']
    allinfo.head()

    fs = 'ENSG00000073792.16_ENSG00000135426.16'
    allinfo['is_fusion']=allinfo['Gene'].apply(
        lambda a: 'IGF2BP2:TESPA1\nfound' if a == fs else 
        'IGF2BP2:TESPA1\nnot found')

    allinfo2 = allinfo
    allinfo = allinfo[allinfo['Gene']==fs]

    a = allinfo.groupby(['Patient','CellType','BarcodeTotal']).sum('IGF2BP2-TESPA1')
    a = a.reset_index()
    b = pd.DataFrame([['Patient 1','None','NoneA',0],['Patient 3','None','NoneB',0]], columns=['Patient','CellType','Barcode','IGF2BP2-TESPA1'])
    a = a.append(b, ignore_index=True)
    a = a.sort_values('Patient')

    plt.figure(figsize=(3, 3), dpi = 600) 
    sns.violinplot(data = a, x = 'Patient', y='IGF2BP2-TESPA1', palette = ['white'])
    sns.swarmplot(data = a, x = 'Patient', y='IGF2BP2-TESPA1', hue = 'CellType', size = 3, 
                  palette = ['grey','#A034F0','#FFA501','#FF2500','#A52A2A','#0433FF'])
    plt.xlabel('', fontsize = 0)
    plt.ylabel('IGF2BP2::TESPA1 expression', fontsize = 11)


    handles =[]
    for i in range(len(labels)):
        handles.append(
            mpatches.Patch(color=colors[i], label=labels[i]))


    plt.legend(handles, labels, ncol=1, frameon=True, title=r'$\bf{Cell\ types}$',loc='center left', bbox_to_anchor=(0.52, 0.78), fontsize = 'xx-small')

    plt.tight_layout()

    plt.savefig(args.fig5c) 

    allinfo = allinfo2
    allinfo = allinfo[allinfo['Patient']=='Patient 2']
    bc_to_isfusion = dict(zip(allinfo['BarcodeTotal'], allinfo['is_fusion']))

    c = Counter(allinfo['BarcodeTotal'])
    df = pd.DataFrame({k:[v] for k,v in c.items()}).T
    df.rename(columns = {0:'n_umi'}, inplace = True)
    df['is_fusion'] = df.index.map(bc_to_isfusion)
    df= df.sort_values('is_fusion')

    plt.figure(figsize=(2.5, 3), dpi = 600) 
    ax = sns.boxplot(data = df, x = 'is_fusion', y='n_umi', palette = ["#FF00FF","#800080"], boxprops=dict(alpha=.5))

    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
        
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
        
        
    plt.xlabel('', fontsize = 0)
    plt.ylabel('Reads per cell', fontsize = 11)

    plt.yscale('log')


    colors = ["#FF00FF","#800080"] 
    labels=['IGF2BP2:TESPA1\nfound','IGF2BP2:TESPA1\nnot found']

    handles =[]
    for i in range(len(labels)):
        handles.append(
            mpatches.Patch(facecolor=colors[i], label=labels[i], alpha = .5, edgecolor='black'))


    plt.legend(handles, labels, ncol=1, frameon=True,loc='center left', bbox_to_anchor=(0.27, 0.88), fontsize = 'x-small')
    plt.tight_layout()
    plt.savefig(args.fig5d) 



def parse_args():
    parser = argparse.ArgumentParser(
        prog='fig5c_fusion_cells.py', 
        usage='python3 fig5c_fusion_cells.py [options]',
        description='violin/swarmplot of cells expressing the IGF2BP2::TESPA1 fusion'
    )
    parser.add_argument(
        '--allinfo', type=str, default = './AllInfo',
        help='AllInfo (scIsoPrep output) file in tsv'
    )
    parser.add_argument(
        '--fig5c', type=str, default='fig5c_violin_plot_fusion.png',
        help='Output png'
    )
    parser.add_argument(
        '--fig5d', type=str, default='fig5d_boxplot_fusion_expression.png',
        help='Output png'
    )


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)

