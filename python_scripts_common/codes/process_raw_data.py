# projects/E-MTAB-8559/python_scripts/codes/process_raw_data.py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread

import matplotlib.pyplot as plt

#### Custom ####

from config import *
from codes.cell_annotation import get_cell_type, annotate_cells, set_obs_colors, add_gene_binary_columns

### Make palettes ### 
#CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)
#GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CLEANED_GENES,markers_file_path=MARKERS_FILE_PATH, adata=None)

#### Functions ####

def import_raw_data_mtx(project_name, raw_data_dict,raw_data_dir,adata_init_path, meta_df_file_path,gene_map_file_path):
    """adata_init  = import_raw_data_mtx(project_name, raw_data_dict,raw_data_dir,adata_init_path, meta_df_file_path,gene_map_file_path=ensembl_to_hgnc_map)"""

    mtx_file = os.path.join(raw_data_dir, raw_data_dict["mtx_filename"])
    cells_file = os.path.join(raw_data_dir, raw_data_dict["cells_filename"])
    genes_file = os.path.join(raw_data_dir, raw_data_dict["genes_filename"])

    # Load matrix
    adata = sc.read_mtx(mtx_file).T
    genes = pd.read_csv(genes_file, sep="\t", header=None)
    cells = pd.read_csv(cells_file, header=None)

    adata.var_names = genes[0].astype(str)
    adata.obs_names = cells[0].astype(str)

    #  Gene mapping 

    gene_map = pd.read_csv(gene_map_file_path)

    if {"ensembl_id", "hgnc_symbol"}.issubset(gene_map.columns):

        ensg_to_hgnc = dict(zip(gene_map["ensembl_id"], gene_map["hgnc_symbol"]))

        adata.var["ensembl_id"] = adata.var_names.str.split(".").str[0]
        adata.var["hgnc_symbol"] = adata.var["ensembl_id"].map(ensg_to_hgnc)
        adata.var_names = adata.var["hgnc_symbol"].fillna(adata.var["ensembl_id"])

    adata.var_names_make_unique()
    adata.var.index.name = None

    #  Metadata 

    meta_df = pd.read_csv(meta_df_file_path)
    adata.obs["barcode"] = adata.obs.index
    adata.obs["donor"] = adata.obs["barcode"].str.split("-").str[0]

    adata.obs = adata.obs.merge(meta_df, on="donor", how="left")

    #  Store raw counts 

    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()

    adata.write(adata_init_path)

    return adata

    
def run_qc(adata, sample_id, adata_qc_path, qc_dir, important_genes):
    """adata_after_qc = run_qc(adata, sample_id, adata_qc_path, qc_dir, important_genes)"""

    sc.settings.autosave = False
    sc.settings.autoshow = False

    # QC metrics 
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "Mt-", "mt-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "rps", "rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    
    sc.pp.calculate_qc_metrics( adata, qc_vars=["mt"], percent_top=None,log1p=False, inplace=True,)

    qc_save_dir = os.path.join(qc_dir, "violin_qc_plots")
    os.makedirs(qc_save_dir, exist_ok=True)
    sc.settings.figdir = qc_save_dir

    sc.pl.violin( adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], save=f"{sample_id}_violin_QC", jitter=0.4,  multi_panel=True,)

    #  Cell QC filtering 
    cells_before = adata.n_obs

    qc_mask = (
        (adata.obs["n_genes_by_counts"] >= adata.obs["n_genes_by_counts"].quantile(0.1))
        & (adata.obs["total_counts"] >= adata.obs["total_counts"].quantile(0.1))
        & (adata.obs["pct_counts_mt"] <= adata.obs["pct_counts_mt"].quantile(0.9))
    )

    print(f"Cells before QC: {cells_before}")
    print(f"Cells after QC: {qc_mask.sum()}")
    print(f"Fraction kept: {qc_mask.mean():.2%}")

    adata = adata[qc_mask].copy()

    #  Preserve raw counts 
    adata.layers["counts"] = adata.X.copy()
    genes_before = adata.var_names.copy()

    sc.pp.filter_genes(adata, min_cells=3)

    removed_genes = genes_before.difference(adata.var_names)
    important_missing = [g for g in important_genes if g in removed_genes]

    if important_missing:

        restore = adata.raw[:, important_missing].copy() if adata.raw is not None else None

        if restore is not None:
            adata = sc.concat([adata, restore], axis=1, merge="same")

        adata.var_names_make_unique()

    #  HVG selection 
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", batch_key="sample_id")
    #  Doublet detection 
    #sc.pp.scrublet(adata, batch_key="sample_id")
    #adata = adata[~adata.obs["predicted_doublet"]].copy()

    #  Normalization 
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.layers["log1p"] = adata.X.copy()

    

    # Force important genes to be HVG
    marker_mask = adata.var_names.isin(important_genes)
    adata.var.loc[marker_mask, "highly_variable"] = True

    #  Scaling + PCA 
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca( adata, n_comps=30, svd_solver="arpack", use_highly_variable=True,)

    #  Graph construction, embedding, neighbors 
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0)

    adata = add_gene_binary_columns(adata, important_genes)

    #  Save 
    adata.write(adata_qc_path)

    return adata


def full_cell_annotation(sample, adata, adata_annotated_path, important_genes, markers_file_path, cell_annotation_file_path, cell_type_label_col="predicted_cell_type"):
    # Cell typing AFTER Leiden, on full gene space
    adata = get_cell_type(sample, adata, markers_file_path, cell_type_label_col)
    adata = annotate_cells(adata, cell_annotation_file_path, cell_type_label_col)
    adata = add_gene_binary_columns(adata, important_genes)

    adata.write(adata_annotated_path)

    return adata


def make_umaps(sample, adata, umap_dir, genes_list,
               gene_color_map, cell_type_colors,
               palette="glasbey",
               max_categories=50):

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import scanpy as sc

    os.makedirs(umap_dir, exist_ok=True)

    # ensure categorical obs
    for col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

    # create gene binary columns
    adata = add_gene_binary_columns(adata, genes_list)

    # compute gene positive counts
    gene_counts = {}

    for g in genes_list:
        col = f"{g}_binary"

        if col not in adata.obs:
            continue

        vals = adata.obs[col].astype(str)

        n_pos = (vals == f"{g}_pos").sum()
        n_total = len(vals)

        gene_counts[g] = (n_pos, n_total)

    # iterate obs columns
    for col in adata.obs.columns:

        # skip unwanted columns
        if col == "barcode":
            continue

        if col.startswith("ulm_score"):
            continue

        n_categories = adata.obs[col].nunique()

        if n_categories > max_categories:
            print(f"Skipping {col} ({n_categories} categories)")
            continue

        umap_save_path = os.path.join(umap_dir, sample)
        os.makedirs(umap_save_path, exist_ok=True)

        # assign colors
        set_obs_colors(adata, col, palette, cell_type_colors, gene_color_map)

        legend_loc = "on data" if col == "leiden" else "right margin"

        title = f"{sample}\n{col}"

        if col.endswith("_binary"):

            gene = col.replace("_binary", "")

            if gene in gene_counts:
                n_pos, n_total = gene_counts[gene]
                title = f"{sample}\n{gene} ({n_pos}/{n_total} cells positive)"

        fig, ax = plt.subplots(figsize=(5,5))

        sc.pl.umap(
            adata,
            color=col,
            legend_loc=legend_loc,
            frameon=False,
            ax=ax,
            title=title,
            show=False,
        )

        plt.tight_layout()

        plt.savefig(
            os.path.join(umap_save_path, f"{col}_umap_{sample}.png"),
            dpi=300,
            bbox_inches="tight"
        )

        plt.close(fig)