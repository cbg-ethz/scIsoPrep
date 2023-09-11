import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re
import math
import argparse
import numpy as np
from scipy.stats import mannwhitneyu,ttest_ind

def mirna(mir_file):
    mir29abc_targets = pd.read_csv(mir_file, sep='\t', index_col = 1)
    mir29abc_targets = mir29abc_targets.sort_values('Gene Symbol').drop_duplicates('Gene Symbol', keep='last')
    mir29abc_targets = list(mir29abc_targets['Gene Symbol'])
    return mir29abc_targets

def label_point(x, y, val,z, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val, 'z':z}, axis=1)
    for i, point in a.iterrows():
        if point['z'] == 'No':
            if  point['val'] == 'COL5A2':
                ax.text(point['x']-.15, point['y']+1.5, str(point['val']), size = '7')
            elif  point['val'] == 'FKBP5':
                ax.text(point['x']-.09, point['y']+.7, str(point['val']), size = '7')
            elif abs(point['x'])>0.4 or point['y']>10:
                ax.text(point['x']-.09, point['y']+1.5, str(point['val']), size = '7')

def main(args):
    mir29abc_targets = mirna(args.mirna)

    dapars = pd.read_csv(args.dapars, sep ='\t')

    dapars.rename(columns = {'DeltaPDUI':'Fraction Change'}, inplace = True)

    all_elong = list(dapars[dapars['Reject']=='No'][dapars[dapars['Reject']=='No']['Fraction Change']>0]['Gene'])
    
    dapars['p-value'].replace(0,1e-199,inplace=True)
    dapars['p-value'] =  dapars['p-value'].apply(lambda a: -math.log10(a) if 0<float(a)<1 else 0)
    dapars['miR29abc'] = dapars['Gene'].apply(lambda a: True if a in mir29abc_targets else False)

    plt.figure(figsize=(3, 3), dpi = 600)
    ax = sns.scatterplot(data=dapars[dapars['Reject']=='Yes'], x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'red'}, alpha = 0.2)

    label_point(dapars['Fraction Change'],dapars['p-value'],dapars['Gene'],dapars['Reject'], plt.gca())

    dapars = dapars[dapars['Reject']=='No'] 
    mir29_elong = dapars[dapars['miR29abc']]['Gene']
    other_elong = [i for i in all_elong if i not in list(mir29_elong)]

    sns.scatterplot(data=dapars[dapars['Fraction Change']<0], x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'#2559A7'}, alpha = 1, ax=ax)
    sns.scatterplot(data=dapars[dapars['Fraction Change']>0], x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'#E8998D'}, alpha = 1, ax=ax)

    sns.scatterplot(data=dapars[dapars['miR29abc']],  x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'#6C9A8B'}, alpha = 1, ax=ax)


    line1 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor='#2559A7')
    line2 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor="lightgrey")
    line3 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor='#E8998D')
    line4 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor="#6C9A8B")

    plt.legend((line1, line2, line3, line4), ('Shortened','Non significant','Lengthened', 'miR29+'), numpoints=1,
               loc='upper left', prop={'size': 8})


    ax.set(xlim=(-0.70, 0.70))

    ax.set(ylim=(0, 72))
    ax.set_ylabel('-log10(p-adjusted)')

    ax.axhline(1.055913, ls='--', color='red',  alpha = 0.5, linewidth =.8)
    ax.axvline(0.1, ls='--', color='red',  alpha = 0.5, linewidth =.8)
    ax.axvline(-0.1, ls='--', color='red', alpha = 0.5, linewidth =.8)

    fig = ax.get_figure()
    fig.savefig(args.fig3b, bbox_inches='tight')

    df_div = pd.read_csv(args.norm_expr)

    all_shorten = list(dapars[dapars['Reject']=='No'][dapars[dapars['Reject']=='No']['Fraction Change']<0]['Gene']) +['Total']
    mir29_elong = list(mir29_elong)+['Total']
    other_elong = list(other_elong)+['Total']

    df_mirna = df_div[df_div['gene'].isin(mir29_elong)]
    df_others = df_div[df_div['gene'].isin(other_elong)]
    df_shorten = df_div[df_div['gene'].isin(all_shorten)]

    df_div_mirna = (df_mirna[df_mirna.columns[1:]]+1).div(df_mirna[df_mirna.columns[1:]].iloc[-1])
    df_div_others = (df_others[df_others.columns[1:]]+1).div(df_others[df_others.columns[1:]].iloc[-1])
    df_div_shorten = (df_shorten[df_shorten.columns[1:]]+1).div(df_shorten[df_shorten.columns[1:]].iloc[-1])

    df_div_mirna['gene'] = df_div['gene']
    df_div_mirna = df_div_mirna[df_div_mirna['gene']!='Total']

    df_div_others['gene'] = df_div['gene']
    df_div_others = df_div_others[df_div_others['gene']!='Total']

    df_div_shorten['gene'] = df_div['gene']
    df_div_shorten = df_div_shorten[df_div_shorten['gene']!='Total']

    df_div_mirna['LogFC'] = np.log(df_div_mirna['Mesothelial.cells.TME']/df_div_mirna['Mesothelial.cells.Distal'])
    df_div_others['LogFC'] = np.log(df_div_others['Mesothelial.cells.TME']/df_div_others['Mesothelial.cells.Distal'])
    df_div_shorten['LogFC'] = np.log(df_div_shorten['Mesothelial.cells.TME']/df_div_shorten['Mesothelial.cells.Distal'])

    dic_plot = {}

    dic_plot['LogFC'] = list(df_div_mirna['LogFC']) + list(df_div_others['LogFC']) + list(df_div_shorten['LogFC'])
    dic_plot['3UTR change'] = ['miR-29a +']*len(df_div_mirna['LogFC']) + ['miR-29a -']*len(df_div_others['LogFC']) + ['Shortened']*len(df_div_shorten['LogFC'])

    df_plot = pd.DataFrame(dic_plot)

    plt.figure(figsize=(4, 3), dpi = 600)
    ax = sns.boxplot(data=df_plot, y ='LogFC', x = '3UTR change',
                     palette = {'miR-29a +':'#6C9A8B','miR-29a -':'#E8998D', 'Shortened':'#2559A7'}, zorder=0)
    sns.stripplot(data = df_plot,y ='LogFC', x = '3UTR change',
                  palette = {'miR-29a +':'#6C9A8B','miR-29a -':'#E8998D', 'Shortened':'#2559A7'}, 
                  dodge=False,jitter = 0.2, alpha =.1, zorder=0, ax = ax)
    ax.set_ylabel('Log(FC) gene expr.')
    ax.set_xlabel('')

    U1, p1 = ttest_ind(list(df_div_mirna['LogFC']),list(df_div_others['LogFC']))
    U2, p2 = ttest_ind(list(df_div_mirna['LogFC']),list(df_div_shorten['LogFC']))


    x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    y, h, col = df_plot['LogFC'].max() + 0.2, 0.2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, "P = {:f}".format(p1), ha='center', va='bottom', color=col)


    x1, x2 = 0, 2   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    y, h1, h2, col = df_plot['LogFC'].max() + 0.2, 0.5, 0.7, 'k'
    plt.plot([x1, x1, x2, x2], [y+h1, y+h2, y+h2, y+h1], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h2, "P = {:f}".format(p2), ha='center', va='bottom', color=col)

    plt.xticks([0, 1, 2], ['Elong.\nmiR-29a +', 'Elong.\nmiR-29a -', 'Shortened'])
    ax.set(ylim=(-1.1, 4))
    fig = ax.get_figure()
    fig.savefig(args.fig3d, bbox_inches='tight')



def parse_args():
    parser = argparse.ArgumentParser(
        prog='fig3b_dapars.py', 
        usage='python3 fig3bd_dapars.py [options]',
        description='plot dapars2_LR output and lengthened/shortened expression boxplot'
    )
    parser.add_argument(
        '--mirna', type=str, default = 'Suppl_Table_3.tsv',
        help='mirna targets file (Supplementary Table 3) in tsv'
    )
    parser.add_argument(
        '--dapars', type=str, default = 'Dapars2_LR_output.txt',
        help='Dapars2_LR output in tsv'
    )
    parser.add_argument(
        '--norm_expr', type=str, default='celltype_expr.csv',
        help='normalized expression of each gene per celltype'
    )    
    parser.add_argument(
        '--fig3b', type=str, default='Fig3b_dapars_3utr.png',
        help='Output png'
    )
    parser.add_argument(
        '--fig3d', type=str, default='fig3d_boxplot_mir29abc.png',
        help='Output png'
    )

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
