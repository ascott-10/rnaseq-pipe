# projects/E-MTAB-8559/python_scripts/main_de_analysis.py


#### Libraries ###

import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc
import gseapy as gp

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread

import matplotlib.pyplot as plt

#### Custom ####

from config import *

from codes.de_expression import rank_genes, plot_ranked_genes, correlate_with_gene, run_cluster_signature
from codes.correlation import get_correlations, make_gene_expression_umaps 
from codes.cell_annotation import make_universal_gene_list, init_color_maps, get_cell_type

### Make palettes ### 
CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)
GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CLEANED_GENES,markers_file_path=MARKERS_FILE_PATH, adata=None)


### Process raw data ###




for project in WHOLE_SAMPLE:

    adata_init_path = os.path.join(WORKING_ADATA_DIR, f"{project}_init.h5ad")
    adata_ranked_path = os.path.join(WORKING_ADATA_DIR, f"{project}_ranked.h5ad")
    adata_qc_path = os.path.join(WORKING_ADATA_DIR, f"{project}_after_qc.h5ad")



    if not os.path.exists(adata_ranked_path):
        adata_ranked = rank_genes(project, adata_init_path, adata_ranked_path, gene="CTCFL")
        adata_ranked.write(adata_ranked_path)
    else:
        adata_ranked = sc.read_h5ad(adata_ranked_path)

    #plot_ranked_genes(adata_ranked_path, project, fig_dir = FIGURES_DIR, tables_dir = TABLES_DIR, group_key="CTCFL_binary", group_name="pos", groupby="label", fdr=0.01, log_fold_change=1)
    #correlate_with_gene(adata_init_path, project, fig_dir = FIGURES_DIR, tables_dir = TABLES_DIR, n_top=50, gene="CTCFL",  method="spearman")
    
    adata_init = sc.read_h5ad(adata_init_path)
    adata_after_qc = sc.read_h5ad(adata_qc_path)

    #top_up, top_down = get_correlations(project)

    make_gene_expression_umaps(project, adata_after_qc, SIG_LIST)
    #make_gene_expression_umaps(project, adata_after_qc, top_down)

    cluster_id_ls =["pos", "neg"]

    if len(cluster_id_ls) == 0:
        continue

    for cluster_id in cluster_id_ls:
        res = run_cluster_signature(
            bdata=adata_after_qc,
            cluster_id=cluster_id,
            sample=project,
            outdir=f"results/signatures"
    )

    
        
    for sample in SAMPLE_LIST:

        bdata_init_path = os.path.join(WORKING_ADATA_DIR, f"{sample}_init.h5ad")
        bdata_ranked_path = os.path.join(WORKING_ADATA_DIR, f"{sample}_ranked.h5ad")
        bdata_qc_path = os.path.join(WORKING_ADATA_DIR, f"{sample}_after_qc.h5ad")

        if not os.path.exists(bdata_ranked_path):

            if not os.path.exists(bdata_init_path):
                bdata_init = adata_init[adata_init.obs["sample"] == sample].copy()
                bdata_init.write(bdata_init_path)
            else:
                bdata_init = sc.read_h5ad(bdata_init_path)


        if not os.path.exists(bdata_ranked_path):
            bdata_ranked = rank_genes(sample, bdata_init_path, bdata_ranked_path, gene="CTCFL")
            bdata_ranked.write(bdata_ranked_path)
        else:
            bdata_ranked = sc.read_h5ad(bdata_ranked_path)


        #plot_ranked_genes(bdata_ranked_path, sample, fig_dir = FIGURES_DIR, tables_dir = TABLES_DIR, group_key="CTCFL_binary", group_name="pos", groupby="label", fdr=0.01, log_fold_change=1)
        #correlate_with_gene(bdata_init_path, sample, fig_dir = FIGURES_DIR, tables_dir = TABLES_DIR, n_top=50, gene="CTCFL",  method="spearman")
        
        bdata_after_qc = sc.read_h5ad(bdata_qc_path)
        #top_up, top_down = get_correlations(sample)

        make_gene_expression_umaps(sample, bdata_after_qc, SIG_LIST)
        #make_gene_expression_umaps(sample, bdata_after_qc, top_down)


        cluster_id_ls =["pos", "neg"]

        if len(cluster_id_ls) == 0:
            continue

        for cluster_id in cluster_id_ls:
            res = run_cluster_signature(
                bdata=bdata_after_qc,
                cluster_id=cluster_id,
                sample=sample,
                outdir=f"results/signatures"
        )

