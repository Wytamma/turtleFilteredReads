#!/usr/bin/bash

# add channels
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# create conda env
conda create --yes --name turtleFilteredReads python=3.6

#install dependancies
echo ''
echo '**************************************************'
echo 'Installing dependancies, this will take a while...'
echo '**************************************************'
echo ''
sleep 3

conda install --yes -n turtleFilteredReads biopython=1.70
conda install --yes -n turtleFilteredReads snakemake=5.1.4

echo ''
echo '**************************************************'
echo 'To activate the virusMAP environment, use'
echo ' $ conda activate virusMAP' 
echo '**************************************************'
echo ''