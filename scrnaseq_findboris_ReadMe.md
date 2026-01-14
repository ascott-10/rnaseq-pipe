# Purpose

This repository contains a modular pipeline for processing and analyzing single-cell RNA-seq (scRNA-seq) data.  
It is designed for use on **NIH Biowulf** with a conda environment but can be adapted to other HPC systems.  


---

# Setup

## 1. Conda environment

```bash
    source myconda
    conda create --name rnaseq-pipe python=3.10
    conda update -n base -c conda-forge conda
    conda activate rnaseq-pipe

    # configure channels
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    # install core dependencies
    conda install pandas numpy matplotlib seaborn scanpy fast_matrix_market anndata
    conda install scvi-tools -c conda-forge         # anndata/scvi-tools support
    conda install -c conda-forge python-igraph      # igraph support
    conda install conda-forge::scanpy 
    conda install conda-forge::decoupler-py         # must be v2.1 or higher
    conda install conda-forge::mudata
    conda install -c conda-forge colorcet

    

*** Setup the folder directories ***
    touch main.py __init__.py config.py .gitignore
    mkdir codes raw_data working_adata results results/figures results/tables
    touch codes/__init__.py codes/utils.py

*** Setup github version-control ***

    git config --global user.name "ascott-10"
    git config --global user.email "ascott10919@gmail.com"

    git init
    git branch -m main
    git add .
    git commit -m "Initial commit"
    git remote add origin <your_repo_url>
    git push -u origin main --force

```
Add a .gitignore to exclude large files, results, and conda environments.

## 2. config.py setup

The pipeline is designed to support the simulataneous analysis of multiple projects at once. The user should assign:
- BASE_DIR
- REFERENCE_DIR
- MARKERS_FILE_PATH
- PROJECT_FOLDERS

in config.py before running `python main.py setup <project_name>`

Analyses include preprocessing raw data, pseudobulk differential expression, transcription factor/pathway scoring,  
and visualization of CTCFL expression.

## 3. Raw data

User should put raw files in the appropriate raw file folder. The pipeline can analyze raw files in many formats including:
- .h5
- .mtx when it includes barcodes, genes, matrix.mtx
- counts in ".csv", ".tsv", ".txt"
- compressed file formats can be unzipped and handled as well

The raw files folder initially will contain the raw files, and from here the sample names will be extracted. Replicates can pose some challenges and the user may have to manually modify the names of replicates to get the exact names they want. 

The raw files folder will also contain the initial, unmodified adata table that is created from the raw files before any processing.


