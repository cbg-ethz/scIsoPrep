#!/usr/bin/env bash

set -Eeuxo pipefail

PWD=$(pwd)

conda create -n scIsoPrep mamba

sed -i "s:<path/to>/scIsoPrep:$PWD:" ./config/*.yaml 

cd scripts
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.1.1.tar.gz
tar -xvf v5.1.1.tar.gz

mamba env update -n scIsoPrep --file SQANTI3.conda_env.yml --prune
conda activate scIsoPrep

cd SQANTI3-5.1.1/

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P ./utilities/
chmod +x ./utilities/gtfToGenePred 

git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install

mamba install -c conda-forge snakemake

cd ../../../

snakemake -s snake/scisoprep.snake --configfile config/config_retina.yaml --cores 100 --use-conda -pr --conda-create-envs-only

