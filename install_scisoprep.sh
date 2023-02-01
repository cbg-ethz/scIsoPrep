#!/usr/bin/env bash

PWD=$(pwd)

sed -i "s:<path/to>/scIsoPrep:$PWD:" ./config/*.yaml 

cd scripts
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v1.6.tar.gz
tar -xvf v1.6.tar.gz

cd SQANTI3-1.6

conda env create -f SQANTI3.conda_env.yml -n scIsoPrep

conda activate scIsoPrep

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P ./utilities/
chmod +x ./utilities/gtfToGenePred 

git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install

conda install -c conda-forge mamba-1.0.0
mamba install -c bioconda snakemake





