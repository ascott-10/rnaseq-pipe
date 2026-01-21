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
from codes.cell_annotation import make_universal_gene_list, init_color_maps, get_cell_type, make_color_map

### Make palettes ### 
CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)
GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CLEANED_GENES,markers_file_path=MARKERS_FILE_PATH, adata=None)

#### Functions ####

def import_raw_data(sample, adata_init_path, raw_file_type = RAW_FILE_TYPE, raw_data_dir = RAW_DATA_DIR, meta_file_path = None):
    
    if raw_file_type == ".h5ad":


        adata_init = sc.read_h5ad(os.path.join(raw_data_dir, f"{sample}.h5ad"))

        adata_init.write(adata_init_path)

        return adata_init
    
    elif raw_file_type == ".mtx":

        mtx_file = glob.glob(os.path.join(RAW_DATA_DIR, "*matrix.mtx"))[0]
        genes_file = glob.glob(os.path.join(RAW_DATA_DIR, "*genes.tsv"))[0]
        cells_file = glob.glob(os.path.join(RAW_DATA_DIR, "*barcodes.tsv"))[0]

        # Load matrix
        adata = sc.read_mtx(mtx_file).T

        # Load genes (23,284 rows)
        genes = pd.read_csv(
            genes_file,
            sep="\t",
            header=None
        )

        # Load cells (19,880 rows)
        cells = pd.read_csv(cells_file, header=None)

        # Assign correctly:
        adata.var_names = genes[0].astype(str)   # genes to var
        adata.obs_names = cells[0].astype(str)   # cells to obs

        print("Genes assigned:", adata.var_names[:5])
        print("Cells assigned:", adata.obs_names[:5])



        #Fix gene names
        gene_map = pd.read_csv(ALL_GENES_FILE_PATH)  # must have: ensembl_gene_id, hgnc_symbol

        if "hgnc_symbol" in gene_map:

            ensg_to_hgnc = dict(zip(gene_map["ensembl_gene_id"], gene_map["hgnc_symbol"]))
            adata.var["ensembl_id"] = adata.var_names.astype(str)
            adata.var["gene_symbol"] = adata.var["ensembl_id"].map(ensg_to_hgnc)
            adata.var["gene_symbol"] = adata.var["gene_symbol"].fillna(adata.var["ensembl_id"])
            adata.var_names = adata.var["gene_symbol"].astype(str)

        adata.var_names_make_unique()
        adata.var.index.name = None
        print(adata.var_names[:10])

        if meta_file_path is not None:
            meta = pd.read_csv(meta_file_path, index_col=0)
            print(meta)

            adata.obs["obs_short"] = [s.split("-")[0] for s in adata.obs_names]

            meta_subset = meta.loc[adata.obs["obs_short"]]

            meta_subset.index = adata.obs.index
            adata.obs = adata.obs.join(meta_subset)


            for col in meta.columns:
                adata.obs[col] = adata.obs[col].astype(str).fillna("")

            adata.obs["sample"] = adata.obs["donor_id"]
            print(adata.obs["sample"].value_counts())

        else:
            adata.obs["sample"] = str(sample)

        adata.obs = adata.obs.apply(lambda s: s.astype(str))

        adata.raw = adata.copy()

        adata.write(adata_init_path)
        return adata




    
    
def run_qc(sample, adata_init_path, adata_qc_path, important_genes = SIG_LIST):

    sc.settings.figdir = FIGURES_DIR
    sc.settings.autosave = False
    sc.settings.autoshow = False

    
    adata = sc.read_h5ad(adata_init_path)

    adata.var['mt'] = adata.var_names.str.startswith(('Mt-', 'MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "rps", "rpl"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    qc_save_dir = os.path.join(FIGURES_DIR, "violin_qc_plots")
    os.makedirs(qc_save_dir, exist_ok=True)

    sc.settings.figdir = qc_save_dir

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        save=f"{sample}_QC",
        jitter=0.4,
        multi_panel=True
    )

    # QC
    cells_before = adata.n_obs

    qc_mask = (
        (adata.obs["n_genes_by_counts"] >= adata.obs["n_genes_by_counts"].quantile(0.1)) &
        (adata.obs["total_counts"]      >= adata.obs["total_counts"].quantile(0.1)) &
        (adata.obs["pct_counts_mt"]     <= adata.obs["pct_counts_mt"].quantile(0.9))
    )

    print(f"Cells before QC: {cells_before}")
    print(f"Cells after QC: {qc_mask.sum()}")
    print(f"Fraction kept: {qc_mask.mean():.2%}")

    adata = adata[qc_mask].copy()

    sc.pp.filter_genes(adata, min_cells=3)

    sc.pp.scrublet(adata, batch_key="sample")
    adata = adata[~adata.obs["predicted_doublet"]].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVGs (flag only)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat_v3"
    )

    # Force important genes to be treated as HVG
    marker_mask = adata.var_names.isin(important_genes)
    adata.var.loc[marker_mask, "highly_variable"] = True

   
    # Scale + PCA using HVGs only
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack", use_highly_variable=True)

    # Neighbors + UMAP + Leiden
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0)

    # Cell typing AFTER Leiden, on full gene space
    adata = get_cell_type(adata)

    if adata.raw is None:
        raise ValueError("adata.raw is None")

    X_raw = adata.raw.X
    genes = adata.raw.var_names.astype(str)
    cells = adata.raw.obs_names.astype(str)

    if sp.issparse(X_raw):
        X_raw = X_raw.toarray()

    df_raw = pd.DataFrame(X_raw, index=cells, columns=genes)

    threshold = 0.1

    for g in ["CTCF", "CTCFL"]:
        if g not in df_raw.columns:
            continue

        expr = df_raw[g].values
        pos_mask = expr > threshold

        adata.obs[f"{g}_binary"] = pd.Categorical(
            np.where(pos_mask, "pos", "neg"),
            categories=["neg", "pos"]
        )

    adata.write(adata_qc_path)

    return adata


def make_umaps(sample, adata_working_path, umap_dir = UMAP_DIR, important_genes=SIG_LIST):
    
    bdata = sc.read_h5ad(adata_working_path)
    if bdata.n_obs == 0:
        return

    
    os.makedirs(umap_dir, exist_ok=True)

    cats_cell_type = bdata.obs["cell_type"].astype(str).unique().tolist()
    cats_leiden = bdata.obs["leiden"].astype(str).unique().tolist()

    cell_type_colors = make_color_map(cats_cell_type)
    leiden_colors = make_color_map(cats_leiden)

    palette_cell_type = {c: cell_type_colors.get(c, "#d3d3d3") for c in cats_cell_type}
    palette_leiden = {c: leiden_colors.get(c, "#d3d3d3") for c in cats_leiden}

    # Use raw counts for full gene list
    X_raw = bdata.raw.X
    genes = bdata.raw.var_names.astype(str)
    cells = bdata.obs_names.astype(str)

    if sp.issparse(X_raw):
        X_raw = X_raw.toarray()

    df = pd.DataFrame(X_raw, index=cells, columns=genes)


    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(
        bdata,
        color="leiden",
        legend_loc="on data",
        frameon=False,
        palette=palette_leiden,
        ax=ax,
        title=f"{sample}\nLeiden Clusters",
        show=False,
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(umap_dir, f"{sample}_leiden.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    sc.pl.umap(
        bdata,
        color="cell_type",
        legend_loc="right margin",
        frameon=False,
        palette=palette_cell_type,
        ax=ax,
        title=f"{sample}\nPredicted Cell Type",
        show=False,
    )
    
    plt.savefig(
        os.path.join(umap_dir, f"{sample}_cell_type.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    genes = [g for g in important_genes if g in bdata.raw.var_names]
    threshold = 0.1

    for g in genes:
        expr = bdata.raw[:, g].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()

        pos_mask = expr > threshold

        bdata.obs[f"{g}_binary"] = pd.Categorical(
            np.where(pos_mask, f"{g}_pos", f"{g}_neg"),
            categories=[f"{g}_neg", f"{g}_pos"],
        )

        expr = df[g].values

        
        threshold = 0.1  
        pos_mask = expr > threshold

        n_pos, n_total = int(pos_mask.sum()), expr.shape[0]

        if g not in GENE_COLOR_MAP:
            new_color = make_color_map([g])[g]
            GENE_COLOR_MAP[g] = new_color

        palette = {
            f"{g}_neg": "#d3d3d3",
            f"{g}_pos": GENE_COLOR_MAP[g],
        }

        fig, ax = plt.subplots(figsize=(5, 5))
        sc.pl.umap(
            bdata,
            color=f"{g}_binary",
            frameon=False,
            palette=palette,
            ax=ax,
            title=f"{sample}\n{g} ({n_pos}/{n_total} cells positive)",
            show=False,
        )
        plt.tight_layout()
        plt.savefig(
            os.path.join(umap_dir, f"{sample}_umap_{g}_binary.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)
