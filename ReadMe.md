# Purpose

This collection of scripts acts as a pipeline for dealing with sc-RNAseq data. It assumes the user is working within the NIH Biowulf and conda environment.

The very first steps involve initiating conda, creating a new environment and updating conda, which may or may not apply to all users. Next conda should be configured with bioconda and conda-forge. FInally, the user should install the basic packages as listed below from the conda environment

# Setup

*** Setup the conda environment ***
    source myconda
    conda create --name rnaseq-pipe python=3.10
    conda update -n base -c conda-forge conda
    conda activate rnaseq-pipe

    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    conda install pandas numpy matplotlib seaborn scanpy fast_matrix_market anndata
    conda install scvi-tools -c conda-forge (for anndata/scvi-tools) 

*** Setup the folder directories ***
    touch main.py __init__.py config.py
    mkdir codes raw_data working_adata results results/figures results/tables
    touch codes/__init__.py

*** Setup github version-control ***

    git config --global user.name "ascott-10"
    git config --global user.email "ascott10919@gmail.com"
    git init
    git branch -m main

# 