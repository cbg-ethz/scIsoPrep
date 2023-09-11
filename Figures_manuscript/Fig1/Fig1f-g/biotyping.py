import pandas as pd
import argparse
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.patches as mpatches
import seaborn as sns


def get_bed_gencode(bed_file):
    bed = pd.read_csv(bed_file, names = list(range(9)), header=None, skiprows = 5, sep = '\t')
    bed.columns = bed.columns.astype(str)
    bed = bed[bed['2'] == 'transcript']
    bed['enst'] = [i.split('transcript_id')[1].split('"')[1] for i in list(bed['8'])]
    bed['biotype'] = [i.split('transcript_type')[1].split('"')[1] for i in list(bed['8'])]

    return dict(zip(bed['enst'],bed['biotype']))

def get_novel_biotypes(novel_biotype):
    bt = pd.read_csv(novel_biotype)
    return dict(zip(bt['pb'],bt['transcript_biotype']))

def legend(Ctypes,CellsColors2,ax):
    handles = []
    labels = []
    for i in range(len(Ctypes)):
        labels.append(Ctypes[i])
        handles.append(
            mpatches.Patch(color=CellsColors2[i], label=Ctypes[i]))
        
    
    ax.legend(handles, labels, ncol=1, frameon=True, title=r'$\bf{Biotype}$',loc='center left', bbox_to_anchor=(1, 0.5))
        
def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{:.1f}%({v:,})'.format(pct,v=val)
    return my_format


def main(args):
    Labels = ["lncRNA":"lncRNA",
     "retained_intron":,"Retained Intron",
      "protein_coding":'Protein Coding',
       "nonsense_mediated_decay":"Nonsense Mediated Decay"]

    biotypes = get_bed_gencode(args.genbed)
    novel_biotypes = get_novel_biotypes(args.novel)

    dict_plot = {'Others':0}
    c = Counter(biotypes.values)
    for bt in c:
        if bt in Labels.keys():
            dict_plot[Labels[bt]] = c[bt]
        else:
            dict_plot["Others"] += c[bt]      

    Colors = sns.color_palette("RdPu", 5, desat=.9)[:5]



    ###Plot was manually done based on v

    fig1, ax1 = plt.subplots(1,2)
    Data = [10095,8036,18560, 7050]
    Values = ['23.1%\n(10,095)','18.4%\n(8,036)','42.4%\n(18,560)','16.1%\n(7,050)']
    ax1[0].pie(Data,
              colors = Colors,
                labels = Values,
               labeldistance=1.2)

    Data_ref   = [75759,29657,85269,17378,24054]
    Values = ['32.6%\n(75,759)','12.8%\n(29,657)','36.7%\n(85,269)','7.5%\n(17,378)','10.4%\n(24,054)']
    ax1[1].pie(Data_ref,
              colors = Colors,
                labels = Values,
              labeldistance=1.2)

    legend(Labels,Colors,fig1)



    fig1.tight_layout(h_pad=-4, w_pad=3)


    plt.savefig(args.out,
                dpi = 600, bbox_inches='tight') 


def parse_args():
    parser = argparse.ArgumentParser(
        prog='fig1_f-g.py', 
        usage='python3 fig1_f-g.py  -o <out.png> [options]',
        description='compute biotypes of gencode db and novel isoforms'
    )
    parser.add_argument(
        '--genbed', type=str, default = '/cluster/work/bewi/members/dondia/projects/ovarian_cancer/reference/gencode.v36.annotation.gtf',
        help='Absolute or relative path(s) to closest to gencode gtf data in bed format'
    )
    parser.add_argument(
        '--novel', type=str, default = './annotation_biotype_analysis_table.csv',
        help='GENCODE defined biotypes of novel isoforms'
    )
    parser.add_argument(
        '-o', '--out', type=str, default='./fig1f-g.png',
        help='Output biotypes of all molecules'
    )


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)



        

