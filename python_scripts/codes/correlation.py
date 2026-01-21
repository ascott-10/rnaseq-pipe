# projects/E-MTAB-8559/python_scripts/codes/correlation.py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread
from scipy.stats import spearmanr

import matplotlib.pyplot as plt

#### Custom ####

from config import *
from codes.cell_annotation import make_universal_gene_list, init_color_maps, make_color_map

### Make palettes ### 
CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)
GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CLEANED_GENES,markers_file_path=MARKERS_FILE_PATH, adata=None)

### Functions ###

def get_correlations(sample, tables_dir=TABLES_DIR, n_genes=10):

    corr_results_dir = os.path.join(tables_dir, "DE_analysis")
    corr_file = os.path.join(
        corr_results_dir, f"{sample}_CTCFL_correlations.csv"
    )

    corr_df = pd.read_csv(corr_file)

    # drop CTCFL itself
    corr_df = corr_df[corr_df["gene"] != "CTCFL"]

    # positive correlations
    corr_pos = (
        corr_df[corr_df["correlation"] > 0]
        .sort_values("correlation", ascending=False)
        .head(n_genes)
    )

    # negative correlations
    corr_neg = (
        corr_df[corr_df["correlation"] < 0]
        .sort_values("correlation", ascending=True)
        .head(n_genes)
    )

    return corr_pos["gene"].tolist(), corr_neg["gene"].tolist()


def make_gene_expression_umaps(sample, adata, gene_list, umap_dir=UMAP_DIR):

    umap_save_path = os.path.join(umap_dir, "gene_expression")
    os.makedirs(umap_save_path, exist_ok=True)

    threshold = 0.1

    for g in gene_list:
        if g not in adata.raw.var_names:
            continue

        expr = adata.raw[:, g].X
        if sp.issparse(expr):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()

        pos_mask = expr > threshold
        n_pos = int(pos_mask.sum())
        n_total = expr.shape[0]

        key = f"{g}_binary"
        adata.obs[key] = pd.Categorical(
            np.where(pos_mask, f"{g}_pos", f"{g}_neg"),
            categories=[f"{g}_neg", f"{g}_pos"],
        )

        if g not in GENE_COLOR_MAP:
            GENE_COLOR_MAP[g] = make_color_map([g])[g]

        palette = {
            f"{g}_neg": "#d3d3d3",
            f"{g}_pos": GENE_COLOR_MAP[g],
        }

        fig, ax = plt.subplots(figsize=(5, 5))
        sc.pl.umap(
            adata,
            color=key,
            frameon=False,
            palette=palette,
            ax=ax,
            title=f"{sample}\n{g} ({n_pos}/{n_total} cells positive)",
            show=False,
        )

        plt.tight_layout()
        plt.savefig(
            os.path.join(umap_save_path, f"{sample}_umap_{g}_binary.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)

        # cleanup
        del adata.obs[key]
