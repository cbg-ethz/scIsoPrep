#!/usr/bin/env bash

set -Eeuxo pipefail

PWD=$(pwd)

sed -i "s:<path/to>/scIsoPrep:$PWD:" ./config/*.yaml 

cd scripts
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.1.1.tar.gz
tar -xvf v5.1.1.tar.gz

cd SQANTI3
sed -i 's/numpy/numpy=1.19.5/' SQANTI3.conda_env.yml 
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

