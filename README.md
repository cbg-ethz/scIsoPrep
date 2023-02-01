# scIsoPrep
A [Snakemake](https://snakemake.github.io/) pipeline for analyzing multiplexed single-cell PacBio concatenated long-reads.

scIsoPrep offers the possibility to unconcatenate, trim, demultiplex large single-cell Pacbio multisample datasets using [IsoSeq3](https://isoseq.how/). It can also collapse transcripts using [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) and classify them using [SQANTI3](https://github.com/ConesaLab/SQANTI3). scIsoPrep first collapses transcripts and filter them per cell, and then repeat this step on all cells together in order to create a common isoforms catalog, using reads attached to isoforms passing all filters in individual cells. 

# Contents
- [Installation](#Installation)
- [Usage](#Usage)
- [Example data](#Example-data)

# Requirements
- Python 3.X
- Conda

# Installation

## Clone repository
First, download scIsoPrep from github and change to the directory:
```bash
git clone https://github.com/cbg-ethz/scisoprep
cd scisoprep
```

## Create conda environment
First, create a new environment named "scIsoPrep" and install all dependencies:
```bash
./install_scisoprep.sh
```

Type `yes` when asked 

# Usage

Before each usage, you should source the scisoprep environment:

```bash
conda activate scIsoPrep
```

The scIsoPrep wrapper script `run_scisoprep.py` can be run with the following shell command:
```bash
./run_scIsoPrep 
```
## Before running the pipeline


* **config file**
  * input directory
    Before running the pipeline, the `config/config.yaml` file needs to be adapted to contain the path to input bam files. It is provided in the first section (`specific`) of the config file.
  * resource information
    In addition to the input path, further resource information must be provided in the section `specific`. This information is primarily specifying
     the genomic reference used for the reads mapping and the transcriptomic reference required for isoform classification. An example `config.yaml` file ready for adaptation, as
    well as a brief description of the relevant config blocks, is provided in the directory `config/`.

* **reference files**
  * A genome fasta file (http://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38)
  * A GENCODE gene annotation gtf file (https://www.gencodegenes.org/human/)

* **sample map**
  * Provide a sample map file, i.e. a tab delimited text file listing all samples that should be analysed, and how many bam files are associated to it (see example below). ID will be used to name files and identify the sample throughout the pipeline.
  * Sample map example:
  ```
  sample     files
  SampleA     2
  SampleB     4
  SampleC     2
  ```
* **input data**
  * This pipeline take as input either concatenated or unconcatenated reads PacBio CCS bam files. I you use concatenated reads input, files should be named `SampleA_1.bam`, `SampleA_2.bam`, `SampleB_1.bam`, etc. (sample name should correspond to the sample map).  If you use unconcatenated reads as input, files should be named `SampleA_1.subreads.bam`, etc.


# Example data

### Downloading data

Data will soon be available on Zenodo

### Run

Edit the `config_retina.yaml` file in the `config` folder, and change the data_path to the folders where bam input files are. Do the same for genome and annotation folder(s). 

Then, edit the `run_scisoprep.py` file and change `--configfile` to `config/config_retina.yaml`.

You can then run the pipeline:
```bash
./run_scIsoPrep 
```

It should run for less thant 2 hours on HPC, and the output file `AllInfo` should be found in the `example_data/results` folder.


