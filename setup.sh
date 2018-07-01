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

mkdir -p tools/
cd tools
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST/ncbi-magicblast-1.3.0-x64-linux.tar.gz -O magicblast.tar.gz
# mac
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST/ncbi-magicblast-1.3.0-x64-macosx.tar.gz -O magicblast.tar.gz

tar  -xvf magicblast.tar.gz
rm -f magicblast.tar.gz
cd ..

echo ''
echo '**************************************************'
echo 'To activate the turtleFilteredReads environment, use'
echo ' $ conda activate turtleFilteredReads' 
echo '**************************************************'
echo ''