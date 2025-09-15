#python_scripts/main.py

#### Libraries ###
import os
import gzip
import shutil

import pandas as pd
import scanpy as sc


import decoupler as dc
import scanpy as sc, matplotlib.pyplot as plt, matplotlib as mpl

sc.set_figure_params(figsize=(6, 6), dpi=150, frameon=False)  # bigger canvas
sc.settings.autoshow = False                                  # don't pop/show
mpl.rcParams['savefig.bbox'] = 'tight'                        # crop safely
mpl.rcParams['savefig.pad_inches'] = 0.1

from config import *
from codes.process_raw_data import detect_multiple_samples, detect_compressed_files, import_raw_data, make_adata
from codes.utils import subset_adata, add_list_to_config
from codes.sc_enrichment import tf_scoring, pathway_scoring, hallmark_genesets, ctcfl_pathway_matrices, ctcfl_geneset_matrices
#### Main ####

def main():
        #1) Detect compressed files and sample names
    detect_compressed_files(raw_data_dir= RAW_DATA_DIR)
    ALL_SAMPLES = detect_multiple_samples(raw_data_dir= RAW_DATA_DIR)
    
    #2) For each sample, import raw data, make AnnData object, preprocess and save working AnnData object
    for sample in ALL_SAMPLES:
        print(f"Processing sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample, raw_file_type = RAW_FILE_TYPE, raw_data_dir = RAW_DATA_DIR)
        adata_working_path = os.path.join(WORKING_DIR, f"{sample}.h5ad")
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir = FIGURES_DIR)
   
        subset_adata(adata, sample, tables_dir=TABLES_DIR)

        #3) For each sample, do TF and pathway scoring, CTCFL pathway and geneset matrices, Hallmark genesets
        adata = tf_scoring(adata, sample, tf_list = ["SMAD1", "SMAD2", "SMAD3", "SMAD4"], plot_yes = True)
        adata = pathway_scoring(adata, sample, path_list = ["NFkB", "TGFb", "Hypoxia"], plot_yes = True)
        adata = ctcfl_pathway_matrices(adata, sample)
        adata = hallmark_genesets(adata, sample, plot_yes = True)
        adata = ctcfl_geneset_matrices(adata, sample)


def working_main():
    ALL_SAMPLES = detect_multiple_samples(raw_data_dir= RAW_DATA_DIR)
    for sample in ALL_SAMPLES:
        print(f"Processing sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample, raw_file_type = RAW_FILE_TYPE, raw_data_dir = RAW_DATA_DIR)
        adata_working_path = os.path.join(WORKING_DIR, f"{sample}.h5ad")
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir = FIGURES_DIR)
        
        adata = tf_scoring(adata, sample, tf_list = ["SMAD1", "SMAD2", "SMAD3", "SMAD4"], plot_yes = True)
        adata = pathway_scoring(adata, sample, path_list = ["NFkB", "TGFb", "Hypoxia"], plot_yes = True)
        adata = ctcfl_pathway_matrices(adata, sample)
        adata = hallmark_genesets(adata, sample, plot_yes = True)
        adata = ctcfl_geneset_matrices(adata, sample)


if __name__ == "__main__":
    ALL_SAMPLES = detect_multiple_samples(raw_data_dir= RAW_DATA_DIR)
    for sample in ALL_SAMPLES:
        print(f"Processing sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample, raw_file_type = RAW_FILE_TYPE, raw_data_dir = RAW_DATA_DIR)
        adata_working_path = os.path.join(WORKING_DIR, f"{sample}.h5ad")
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir = FIGURES_DIR)
        subset_adata(adata, sample, tables_dir=TABLES_DIR)

    
#End of main.py