# Purpose



This repository contains a modular, project-isolated pipeline for processing and analyzing single-cell RNA-seq (scRNA-seq) data.

It is designed for use on NIH Biowulf with conda and Slurm, but can be adapted to other HPC systems.

Key design principles:

* each project lives in its own directory
* project-specific logic lives with the project
* a single shared reference dataset
* reproducible, restartable execution
* Slurm-safe at every step

# Setup

## 1.  Conda environment

This project requires Python 3.10.

Create and activate the environment:

    conda create -n rnaseq-pipe python=3.10
    conda activate rnaseq-pipe

From the root of the repository (where requirements.txt lives):

    pip install -r requirements.txt

## 2. Repository strucutre

    scrnaseq_pub/
    ├── projects/               # Each proect is fully contained within its own project folder
    │   ├── human_fetal_gonad/
    │   │   ├── raw_data/       # Raw data (*.mtx, *.tsv, *.h5ad, *h5) formats; single or multi-sample files
    │   │   ├── working_adata/  # Procesed raw data (*.h5ad)
    │   │   ├── results/        # Figures, tables
    │   │   └── python_scripts/ # Contains editable config, main driver script
    │   │       ├── config.py
    │   │       ├── main_process_raw_data.py
    │   │       └── codes/      # Contains codes that populate the main driver
    │   ├── ovarian_cancer/
    │   ├── tubes/
    │   └── ...
    │
    ├── reference_data/         # For storage of global reference files
    │   ├── ensembl_hgnc.csv
    │   ├── ensembl_genes.csv
    │   ├── markers_df.csv
    │   ├── markers_df_fetal_gonad.csv
    │   └── raw_data_sources.xlsx
    │
    ├── run_scrnaseq_fetal_gonad.sh    # Bash/slurm scripts for each project
    ├── run_scrnaseq_MTAB.sh
    ├── run_scrnaseq_ovarian_cancer.sh
    ├── run_scrnaseq_tubes.sh
    │
    ├── logs/
    ├── errors/
    ├── requirements.txt
    └── scrnaseq_findboris_ReadMe.md



## 3. config,py

Each project has its own config.py.

It defines:

* paths (raw data, output directories)
* sample list
* marker files
* genes of interest that are:
       - protected from filtering  
       - retained after HVG selection  
       - used for binary expression calls  
        - visualized on UMAPs  
* project-specific metadata handling

## 4. raw data

Supported formats:

* MTX (10x / GSM)

    Files must be uncompressed and follow 10x conventions:

    *.matrix.mtx  
    *.barcodes.tsv  
    *.genes.tsv  or  *.features.tsv

    For multi-sample MTX projects, files must share a prefix:

    GSM5599221_Norm2.matrix.mtx  
    GSM5599221_Norm2.barcodes.tsv  
    GSM5599221_Norm2.genes.tsv      


* CSV / TSV counts

Compressed files (.gz, .tar.gz) should be unzipped before running.

# Running the pipeline

Each project is run independently.

Local / interactive

From the project directory:

    cd projects/<PROJECT_NAME>
    python python_scripts/main_process_raw_data.py


TSlurm Execution (recommended)

Each project has its own Slurm wrapper.

Example: run_scrnaseq_tubes.sh

    #!/bin/bash
    #SBATCH --job-name=tubes_scrnaseq
    #SBATCH --output=logs/tubes_scrnaseq_%j.out
    #SBATCH --error=errors/tubes_scrnaseq_%j.err
    #SBATCH --time=24:00:00
    #SBATCH --mem=64G
    #SBATCH --cpus-per-task=4

    source /data/scottaa/conda/etc/profile.d/conda.sh
    conda activate rnaseq-pipe

    cd /data/scottaa/scrnaseq_pub/projects/tubes
    python python_scripts/main_process_raw_data.py


Run with:

    sbatch run_scrnaseq_array.sh

# Pipeline Steps 

Each project runs as an independent Slurm task.

For each sample:

1. Raw data import  
2. QC - cell filtering, mitochondrial filtering, scrublets  
3. Normalization, HVG selection, PCA/Neighbors/UMAP/Leiden clustering  
4. Cell-type annotation  
5. Binary expression UMAPs for important genes  


Intermediate results are written to working_adata/ and reused automatically.