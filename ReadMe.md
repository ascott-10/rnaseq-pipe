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
    conda install -c conda-forge python-igraph (for igraph package)
    conda install conda-forge::scanpy 
    conda install conda-forge::decoupler-py (NOTE THIS --> must be vers 2.1 or higher)
    conda install conda-forge::mudata
    

*** Setup the folder directories ***
    touch main.py __init__.py config.py .gitignore
    mkdir codes raw_data working_adata results results/figures results/tables
    touch codes/__init__.py codes/utils.py

*** Setup github version-control ***

    git config --global user.name "ascott-10"
    git config --global user.email "ascott10919@gmail.com"
    git init # initialize git
    git branch -m main # change main branch name to "main"
    git add . # add files to staging
    git commit -m "Initial commit"
    git push -u origin main --force

    *Setup .gitignore

# Scripts

    main.py will contain the main driver functions
    utils.py will contain modular functions that will be used across multiple codes
    process_data.py will contain the initial functions that process raw data into adata tables and perform initial QC and normalization



# Workflow

The workflow begins with initial data, whether that be in preprocessed .h5 or in .mtx/genes/barcodes and then converts them to anndata tables. The adata tables can then be used for downstream expression analysis.

The config file should be setup with the filepaths to the base folders (project, results, working adata). These can be customized.

Raw files (h5) should be placed inside raw_data, and processed adata tables should be placed inside working_adata.


## Step 1

The user should upload the raw data files (in .h5, .mtx or .csv) into raw_data and update .config to indicate which file type. The user should then set the file paths and directories as appropriate. In the main function, the function detects sample names and if any files are compressed, they are then uncompressed before moving on. 

## Step 2

The main script iterates through the samples and creates adata tables and performs normalization, filtering and finds highly variable genes. It also assigns cell types if there aren't any and creates UMAPs highlighting cell types versus CTCFL expression (based on raw counts). These are printed to console; the figures are saved in results/figures/umap by default. 

## Step 3

For each sample, do TF and pathway scoring, CTCFL pathway and geneset matrices, Hallmark genesets. THe user can customize which transcription factors to focus and which pathways to focus on. The user should pay attention to where the plots will be saved and what their names are but generally everything will be saved under figures

