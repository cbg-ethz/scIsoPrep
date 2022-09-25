###misc_snake.py
import os.path
import sys
import pandas as pd

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq)
    letters = [complement[base] for base in bases] 
    letters = ''.join(letters)
    reverse = letters[::-1]
    return reverse

def getSampleNames():
    sample2id = {}
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    df = pd.read_csv(SAMPLEMAPPING, sep = '\t')
    samples = dict(zip(df['sample'],df['files']))
    for sample,id in samples.items():
        sample2id[sample]=list(range(1,id+1))

    return sample2id

def getBCs(samples):
    BCs = {}
    for sample in samples:
        bc_file = '{}/{}.txt'.format(BC_DIR,sample)
        try:
            open(bc_file, "r")
        except IOError:
            print('{} is missing. A barcodes file should be provided for each sample!'.format(bc_file))
            sys.exit(1)
        df_bc = pd.read_csv(bc_file, sep = '\t')
        bc = list(df_bc.barcodes)
        rev_bc = [reverse_complement(i) for i in bc]
        BCs[sample] = rev_bc

    return BCs

def sample2ids(wildcards):
    return expand('{{sample}}/polyA_trimming/{{sample}}_{id}.fltnc.bam', 
               id = IDs[wildcards.sample])

def sample2bcs1(wildcards):
    return expand(os.path.join("{{sample}}/cells/{barcode}_{{sample}}/sqanti3_qc",
        "x.collapsed_classification.filtered_lite_classification.txt"), 
        barcode = BCs[wildcards.sample])
def sample2bcs2(wildcards):
    return expand(os.path.join("{{sample}}/cells/{barcode}_{{sample}}/mapping",
        "x.collapsed.group.txt"),
        barcode = BCs[wildcards.sample])
def sample2bcs3(wildcards):
    return expand(os.path.join("{{sample}}/cells/{barcode}_{{sample}}/mapping",
        "mapped.fasta.sam"), 
        barcode = BCs[wildcards.sample])
def sample2bcs4(wildcards):
    return expand(os.path.join("{{sample}}/cells/{barcode}_{{sample}}/dedup",
        "dedup.info.csv"), 
        barcode = BCs[wildcards.sample])

