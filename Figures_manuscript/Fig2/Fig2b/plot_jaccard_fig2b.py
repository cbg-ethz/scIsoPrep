import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import jaccard_score
import itertools
import argparse

def prep_df(df, index, vals):
    df = pd.pivot(df, index=index, columns=vals, values = vals)
    df.sort_values(index, inplace = True)
    df = df.fillna(0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(1)
    return df

def compute_jaccard(df_short, df_long):
    my_df = pd.DataFrame()
    i=0

    for s_col,l_col in itertools.product(df_short.columns, df_long.columns):
        j = i%7
        ndf = pd.DataFrame({'short':df_short[s_col], 'long' :df_long[l_col]})
        ndf = ndf.loc[~(ndf==0).all(axis=1)]
        jac = jaccard_score(ndf['short'], ndf['long'])
        my_df.loc[j,s_col] = jac
        i+=1
    my_df.index = my_df.columns
    return my_df

def plot_jaccard(my_df,out):
    plt.figure(figsize=(4, 4), 
               dpi = 600)
    
    labels = my_df.applymap(lambda v: "{:.2f}".format(v) if v > 0.1 else '')
    hm = sns.heatmap(my_df, annot=labels, cmap='RdPu',annot_kws={'fontsize':9}, fmt='')
    

    fig = hm.get_figure()
    fig.tight_layout()
    fig.savefig('fig2b_jaccard_heatmap_{}.png'.format(out)) 


def main(args)
    df_iso = pd.read_csv(args.iso)
    df_long = pd.read_csv(args.long)
    df_short = pd.read_csv(args.short)

    df_long = prep_df(df_long,'Barcodes','Cluster')
    df_short = prep_df(df_short,'Barcodes','Cluster')
    df_iso = prep_df(df_iso,'Barcodes','Cluster')



    short_v_long = compute_jaccard(df_short, df_long)
    short_v_iso = compute_jaccard(df_short, df_iso)
    long_v_iso = compute_jaccard(df_long, df_iso)

    plot_jaccard(short_v_long, 'SR_vs_LR')
    plot_jaccard(short_v_iso, 'SR_vs_iso')
    plot_jaccard(long_v_iso, 'LR_vs_iso')


def parse_args():
    parser = argparse.ArgumentParser(
        prog='plot_jaccard_fig2b.py', 
        usage='python3  [options]',
        description='plotting jaccard distance between cell annotation and manual clustering'
    )
    parser.add_argument(
        '--iso', type=str, default='umapCoord_iso.csv',
        help='path to manually defined pacb isoform level umap clustering'
    )
    parser.add_argument(
        '--long', type=str, default='umapCoord_long.csv',
        help='path to manually defined pacb gene level umap clustering'
    )
    parser.add_argument(
        '--short', type=str, default='umapCoord_short.csv',
        help='path to manually defined illumina umap clustering'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)