import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re
import math
import argparse
import numpy as np
from scipy.stats import mannwhitneyu,ttest_ind

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import re
import math
from scipy import stats

def mirna(mir_file):
    mir29abc_targets = pd.read_csv(mir_file, sep='\t', index_col = 1)
    mir29abc_targets = mir29abc_targets.sort_values('Gene Symbol').drop_duplicates('Gene Symbol', keep='last')
    mir29abc_targets = list(mir29abc_targets['Gene Symbol'])
    return mir29abc_targets

def label_point(x, y, val,z, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val, 'z':z}, axis=1)
    for i, point in a.iterrows():
        if point['z'] == 'No':
            if abs(point['x'])>0.5 or point['y']>30:
                ax.text(point['x']-.09, point['y']+1.5, str(point['val']), size = '7')

def main(args)
    dapars = pd.read_csv(args.dapars, sep ='\t')
    mir29abc_targets = mirna(args.mirna)
    dapars.rename(columns = {'DeltaPDUI':'Fraction Change'}, inplace = True)

    dapars = dapars[~dapars['Gene'].str.contains("MT")]
    dapars = dapars[~dapars['Gene'].str.contains("RP")]
    dapars = dapars[~dapars['Gene'].str.contains("H2")]

    all_elong = list(dapars[dapars['Reject']=='No'][dapars[dapars['Reject']=='No']['Fraction Change']>0]['Gene'])
    
    dapars['p-value'].replace(0,1e-199,inplace=True)
    dapars['p-value'] =  dapars['p-value'].apply(lambda a: -math.log10(a) if 0<float(a)<1 else 0)
    dapars['p-value']=  dapars['p-value'].apply(lambda a: 50 if a > 50 else a)
    dapars['miR29abc'] = dapars['Gene'].apply(lambda a: True if a in mir29abc_targets else False)

    plt.figure(figsize=(3, 3), dpi = 600)
    ax = sns.scatterplot(data=dapars[dapars['Reject']=='Yes'], x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'red'}, alpha = 0.2)

    label_point(dapars['Fraction Change'],dapars['p-value'],dapars['Gene'],dapars['Reject'], plt.gca())

    dapars = dapars[dapars['Reject']=='No'] 

    sns.scatterplot(data=dapars[dapars['Fraction Change']<0], x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'#2559A7'}, alpha = 1, ax=ax)
    sns.scatterplot(data=dapars[dapars['Fraction Change']>0], x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'#E8998D'}, alpha = 1, ax=ax)
    sns.scatterplot(data=dapars[dapars['miR29abc']],  x = 'Fraction Change', y='p-value',
                        hue = 'Reject', palette = {'Yes':'grey','No':'#6C9A8B'}, alpha = 1, ax=ax)

    #'#6C9A8B','miR-29a -':'#FFCAB1', 'Shortened':'#2559A7'


    line1 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor='#2559A7')
    line2 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor="lightgrey")
    line3 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor='#E8998D')
    line4 = Line2D([], [], color="white", marker='o', markersize=7, markerfacecolor="#6C9A8B")

    plt.legend((line1, line2, line3, line4), ('Shortened','No signal','Lengthened', 'miR29+'), numpoints=1,
               loc='upper left', prop={'size': 8}, bbox_to_anchor=(-0.015, 1.02))

    ax.set(ylim=(0, 80))
    ax.set_ylabel('-log10(p-adjusted)')

    ax.axhline(1.055913, ls='--', color='red',  alpha = 0.5, linewidth =.8)
    ax.axvline(0.1, ls='--', color='red',  alpha = 0.5, linewidth =.8)
    ax.axvline(-0.1, ls='--', color='red', alpha = 0.5, linewidth =.8)


    fig = ax.get_figure()
    fig.savefig(args.fig4f, bbox_inches='tight')


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
        '--fig4f', type=str, default='Fig4f_dapars_3utr.png',
        help='Output png'
    )


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
