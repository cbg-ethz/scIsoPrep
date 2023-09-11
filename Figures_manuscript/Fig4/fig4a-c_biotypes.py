import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from collections import defaultdict,Counter

import alluvial #https://github.com/vinsburg/alluvial_diagram.git

seed=13
np.random.seed(seed)


def gene_name_dico(file):
    dico = {}
    enst_to_biotype = {}
    names = ['chr','source','type','start','end', '_', 'strand', '__', 'info']
    gencode = pd.read_csv(file, sep = '\t', names = names, skiprows = 5)
    for i,row in gencode.iterrows():
        if row['type'] == 'gene':
            ensg = row['info'].split(';')[0].split('"')[1]
            gene = row['info'].split(';')[2].split('"')[1]
            dico[ensg] = gene
        elif row['type'] == 'transcript':
            enst = row['info'].split('transcript_id')[1].split('"')[1]#.split('.')[0]
            try:
                btype = row['info'].split('transcript_type')[1].split('"')[1]
            except IndexError:
                continue
            enst_to_biotype[enst] = btype

    return dico, enst_to_biotype, gencode

def order_labels(g,d,  HGSOC_counter, right = False):
    tot = sum([i for dic in HGSOC_counter.values() for i in dic.values()])
    labels = []
    if right:
        for bt in d:
            for bt_order in g:
                labels.append(HGSOC_counter[bt_order][bt])
            labels.append(tot*.03)
    else:
        for bt in g:
            for bt_order in d:
                labels.append(HGSOC_counter[bt][bt_order])
            labels.append(tot*.03)
                
    return list(reversed(labels[:-1]))

def sublabels(xs, tot, right = False):
    v = 0
    xpos = 0.01
    ha = 'left'
    if right:
        xpos = 0.98
        ha = 'right'
    for x in xs:
        v+=x
        if x > tot*.03:
            plt.text(xpos,
                     v-x/2,
                     str(x),
                     fontsize = 20,
                     ha=ha, va='center')

def get_biotypes(df,gene_to_maxpi,isoID)
    for idx, row in df.iterrows():
        gene_to_maxpi[row['ENSG']] = [row['maxDeltaPI_ix1'],row['maxDeltaPI_ix2']]

    biotypes = {'HGSOC':[],'Distal':[]}
    protein_coding =0
    HGSOC_protein_coding = 0
    HGSOC_dict = defaultdict(lambda:[])

    rename = {'nonsense_mediated_decay' : 'NMD',  
              'protein_coding':'Protein-coding', 
              'retained_intron': 'Retained intron',
              'processed_transcript': 'Non protein-coding',
              'TEC':'Non protein-coding',
              'processed_pseudogene':'Non protein-coding',
              'snoRNA':'Non protein-coding',
              'lncRNA':'Non protein-coding',
              'transcribed_processed_pseudogene':'Non protein-coding',
              'transcribed_unprocessed_pseudogene':'Non protein-coding',
              'non_stop_decay':'Non protein-coding',
              'misc_RNA':'Non protein-coding',
             } 


    ENSG_isoID_to_biotype = dict(zip(isoID['ENSG_isoID'],isoID['biotype']))
    pi_df = pd.read_csv(args.iso_count, sep = '\t')

    pi_df['ENSG_isoID'] = pi_df[['Gene','IsoID']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    pi_df['biotype'] = pi_df['ENSG_isoID'].map(ENSG_isoID_to_biotype)
    o=0
    for k,v in gene_to_maxpi.items():
        #for i in [0,1]:
        gene = pi_df[pi_df['Gene'] == k]
        HGSOC = gene[gene['delta']==max(gene['delta'])]['biotype'].to_string(index = False).split('\n')[0]
        Distal = gene[gene['delta']==min(gene['delta'])]['biotype'].to_string(index = False).split('\n')[0]
        Distal = Distal.strip()
        HGSOC = HGSOC.strip()
        o+=1
        if HGSOC in rename.keys():
            HGSOC = rename[HGSOC]
        if Distal in rename.keys():
            Distal = rename[Distal]
        
        HGSOC_dict[HGSOC].append(Distal)
        if HGSOC != Distal:
            biotypes['HGSOC'].append(HGSOC)
            biotypes['Distal'].append(Distal)
            #HGSOC_dict[HGSOC].append(Distal)
            if 'Protein-coding' in [HGSOC, Distal]:
                protein_coding +=1
                if HGSOC == 'Protein-coding':
                    HGSOC_protein_coding+=1              

    HGSOC_counter = {k+' ': Counter(v) for k,v in HGSOC_dict.items()}

    return HGSOC_counter


def main(args)
    dico, enst_to_biotype, gencode = gene_name_dico(args.gencode)

    classif = pd.read_csv(args.classif, sep='\t')
    pb_to_enst = dict(zip(classif['isoform'],classif['associated_transcript']))

    novel_biotype = pd.read_csv(args.new_biotypes)
    novelpb_to_biotype = dict(zip(novel_biotype['pb'],novel_biotype['transcript_biotype']))

    isoID = pd.read_csv(args.isoid, sep = '\t')
    isoID['ensg'] = isoID['Gene'].apply(lambda a: dico[a] if a in dico.keys() else a)
    isoID['novelbiotype'] = isoID['Isoform'].apply(lambda a: novelpb_to_biotype[a] if a in novelpb_to_biotype.keys() else a)
    isoID['enst'] = isoID['novelbiotype'].apply(lambda a: pb_to_enst[a] if a in pb_to_enst.keys() else a)
    isoID['biotype'] = isoID['enst'].apply(lambda a: enst_to_biotype[a] if a in enst_to_biotype.keys() else a)
    isoID['biotype'] = isoID['biotype'].apply(lambda a: 'Unknown' if a == 'novel' else a)
    isoID = isoID.fillna('Unknown')
    #isoID['biotype'] = isoID['biotype'].apply(lambda a: 'Non MANE' if a not in ['MANE'] else a)
    #isoID['biotype'] = isoID['biotype'].apply(lambda a: 'NA' if a not in ['1','2','3','4','5'] else a)

    isoID['ENSG_isoID'] = isoID[['Gene','iso.id']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    gene_isoID_to_biotype = dict(zip(isoID['ENSG_isoID'],isoID['biotype']))
    gene_isoID_to_enst = dict(zip(isoID['ENSG_isoID'],isoID['enst']))

    
    df = pd.read_csv(DIU_results, sep = '\t')

    print(len(df))
    df['FDR'] = df['FDR'].apply(lambda x: np.where(x < 1e-300,1e-300,x))
    df['-log10(FDR)'] = -np.log10(df['FDR'])

    df = df.reindex(df['dPI'].abs().sort_values(ascending = False).index)
    df = df.reset_index()
    df['dPI_index'] = df.index
    df = df.reindex(df['FDR'].abs().sort_values(ascending = True).index)
    df = df.reset_index()
    df['FDR_index'] = df.index
    df['ranksum'] = (df['FDR_index']+1) * (df['dPI_index']+1)
    df = df.sort_values('ranksum')
    df = df.fillna(-1)
    df = df.astype({'maxDeltaPI_ix2':'int'})
    df['ENSG_isoID1'] = df[['Gene','maxDeltaPI_ix1']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    df['ENSG_isoID2'] = df[['Gene','maxDeltaPI_ix2']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    df['maxDeltaPI_ix1_biotype'] = df['ENSG_isoID1'].apply(lambda a: gene_isoID_to_biotype[a] if a in gene_isoID_to_biotype.keys() else a)
    df['maxDeltaPI_ix2_biotype'] = df['ENSG_isoID2'].apply(lambda a: gene_isoID_to_biotype[a] if a in gene_isoID_to_biotype.keys() else 'NA')
    df['ENSG'] = df['Gene']
    df = df.replace({"Gene": dico})

    df['ENST1'] = df['ENSG_isoID1'].apply(lambda a: gene_isoID_to_enst[a] if a in gene_isoID_to_enst.keys() else a)
    df['ENST2'] = df['ENSG_isoID2'].apply(lambda a: gene_isoID_to_enst[a] if a in gene_isoID_to_enst.keys() else a)


    df['Gene'] = df['Gene'].apply(lambda a: '{}-{}'.format(dico[a.split('-')[0]],dico[a.split('-')[1]] ) if '-ENSG' in a else a)
    df = df[~df['Gene'].str.contains("MT")]
    df = df[~df['Gene'].str.contains("RP")]
    df = df[~df['Gene'].str.contains("H2")]

    gene_to_maxpi = {}
    lengths = {}
    lengths['all'] = len(df)
    counters ={}
    figname = {'0.2':'Fig4b_dpi0.2_biotypes.png', '0.2':'Fig4c_dpi0.5_biotypes.png'}

    for dpi in [0.2,0.5]:
        df_dpi = df[(df['FDR']<0.05) & (abs(df['dPI']) > dpi)]
        lengths['dPI = {}'.format(str(dpi))] = len(df_dpi)
        HGSOC_counter = get_biotypes(df_dpi,gene_to_maxpi,isoID)
        input_data = HGSOC_counter

        ax = alluvial.plot(
            input_data,  alpha=.4, color_side=0, figsize=(6,8), rand_seed=seed,
            disp_width=True, wdisp_sep=' ', fontname='Arial', fontsize=50, 
             res=19, h_gap_frac=0.045,  width_in=True)

        fig = ax.get_figure()
        #fig.set_size_inches(5,5)

        ax.text(-.3, ax.get_ylim()[1], 'HGSOC', fontsize=36)
        ax.text(1, ax.get_ylim()[1], 'Distal', fontsize=36)

        tot = sum([i for dic in HGSOC_counter.values() for i in dic.values()])

        g = ['Protein-coding ','Non protein-coding ', 'Unknown ','Retained intron ', 'NMD ' ]
        d = ['Protein-coding','Non protein-coding', 'Retained intron', 'Unknown', 'NMD' ]  
        x1= order_labels(g,d,HGSOC_counter)
        x2 = order_labels(g,d,HGSOC_counter, right = True)
        sublabels(x1, tot)
        sublabels(x2, tot, right = True)
        fig.savefig(figname[str(dpi)],bbox_inches='tight')


    species = (
        "HGSOC vs. Distal"
    )
    weight_counts = {
        'No Signal': np.array([lengths['all']-lengths['dPI = 0.2']]),
        ">20% change": np.array([lengths['dPI = 0.2']-lengths['dPI = 0.5']]),
        '>50% change \n(switch)': np.array([lengths['dPI = 0.5']])
    }
    width = 0.5

    fig, ax = plt.subplots(figsize=(1, 3), dpi = 600)
    bottom = np.zeros(3)
    for boolean, weight_count in weight_counts.items():
        p = ax.bar(species, weight_count, width, label=boolean, bottom=bottom)
        bottom += weight_count
    ax.set_title("Change in\n isoforms usage")
    ax.set_ylabel('Genes')
    ax.legend(loc="lower right", bbox_to_anchor=(1.3, -.5), fontsize = 8)
    fig.savefig('Fig4a_DIU.png',bbox_inches='tight')



def parse_args():
    parser = argparse.ArgumentParser(
        prog='fig4a-c_biotypes.py', 
        usage='python3 fig4a-c_biotypes.py [options]',
        description='fig 4a-c DIU barchart plot and biotypes alluvial plots'
    )
    parser.add_argument(
        '--gencode', type=str, default = 'gencode.v36.annotation.gtf',
        help='gencode gtf'
    )
    parser.add_argument(
        '--classif', type=str, default = 'correct_merge.collapsed_classification.txt',
        help='SQANTI output classification'
    )
    parser.add_argument(
        '--new_biotypes', type=str, default='annotation_biotype_analysis_table.csv',
        help='csv containing biotypes of novel isoforms'
    )    
    parser.add_argument(
        '--isoid', type=str, default='Iso-IsoID.csv',
        help='isoID file from scisorseqr'
    )
    parser.add_argument(
        '--DIU_results', type=str, default='HGSOC-all_Distal-all_25_results.csv',
        help='Differential isoforms usage analysis results from scisorseqr'
    )
    parser.add_argument(
        '--iso_count', type=str, default='HGSOC-all_Distal-all_isocount.tsv',
        help='isoform count and fraction per condition from scisorseqr'
    )

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)

