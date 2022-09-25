# scIsoPrep
A [Snakemake](https://snakemake.github.io/) pipeline for analyzing multiplexed single-cell PacBio concatenated long-reads.

scIsoPrep offers the possibility to unconcatenate, trim, demultiplex large single-cell Pacbio multisample datasets using [IsoSeq3](https://isoseq.how/). It can also collapse transcripts using [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) and classify them using [SQANTI3](https://github.com/ConesaLab/SQANTI3). scIsoPrep first collapses transcripts and filter them per cell, and then repeat this step on all cells together in order to create a common isoforms catalog, using reads attached to isoforms passing all filters in individual cells. 

# Contents
- [Installation](#Installation)
- [Usage](#Usage)
- [Example data](#Example-data)

# Requirements
- Python 3.X
- Snakemake

# Installation

## Clone repository
First, download scIsoPrep from github and change to the directory:
```bash
git clone https://github.com/cbg-ethz/scisoprep
cd scisoprep
```

## Create conda environment
First, create a new environment named "scisoprep":
```bash
conda create --name scisoprep python=3
```

Second, source it:
```bash
conda activate scisoprep
```

## Install requirements

scIsoPrep follows the best practices of the Snakemake workflow manager in providing the software needed to run the pipeline in per-rule conda environments. Those environmnents are specified in the `envs/` directory in yaml files that are named `{rule_name}.yaml`. The easiest way to install and use the software is by running Snakemake with the `--use-conda` parameter. Snakemake will try to find the environments of the yaml files the rules point to, and install them if they are not already available. The directory for installing the conda environments can be specified with the `--conda-prefix` parameter.

1. Make sure `snakemake` is in your PATH.
   Follow the instructions on how to install `snakemake` [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
2. Install all conda environments the workflow needs before running an analysis with

```bash
snakemake --use-conda --conda-create-envs-only --conda-prefix /my/directory/for/conda/envs/ -s workflow/snakefile_basic.smk --configfile config/config.yaml
```

- `--use-conda` instructs snakemake to utilize the `conda:` directive in the rules
- `--conda-create-envs-only` specifies that only the installation of conda environments is triggered, not the analysis of the samples.
- *(optional):* with `--conda-prefix /my/directory/for/conda/envs/` a directory for the installation of the conda environments can be specified.


Now you are ready to run **scIsoPrep**!

# Usage
The scIsoPrep wrapper script `run_scisoprep.py` can be run with the following shell command:
```bash
./run_scIsoPrep 
```
### Before running the pipeline

* **config file**
  * input directory
    Before running the pipeline the `config.yaml` file needs to be adapted to contain the **path to input fastq files** for the intended analysis. It is provided in the first section (`specific`) of the config file.
  * resource information
    In addition to the input path, further resource information must be provided in the section `specific`. This information is primarily specifying
     the genomic reference used for the reads mapping and the transcriptomic reference required for isoform classification. An example `config.yaml` file ready for adaptation, as
    well as a brief description of the relevant config blocks, is provided in the directory `config/`.
* **sample map**
Provide a "sample_map", i.e. a tab delimited text file listing all samples that should be analysed (one row per sample).
The sample map must contain a column with the header `sample` (see example below). This ID will be used to name files and identify the sample throughout the pipeline.
An example file ready for adaptation is provided in the directory `config/`.

Sample map example:
```
sample           files
SAMPLE-1_Tumor     2
SAMPLE-1_Healthy   4
SAMPLE-2_Tumor     2
SAMPLE-2_Healthy   4
```

# Example data

TBD