# python_scripts/codes/process_adata.py

#### Libraries ###
import os
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread
import matplotlib.pyplot as plt
import seaborn as sns
import anndata

import decoupler as dc

import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation

# Plot defaults
sc.set_figure_params(figsize=(3, 3), frameon=False)

# Project configuration (paths, constants, etc.)
from config import *
import gzip
import shutil
from config import *

from codes.utils import subset_adata, get_cell_type, add_list_to_config, get_raw_gene_counts, init_color_maps, make_color_map, GENE_COLOR_MAP



init_color_maps(CTCFL_MARKERS + ALL_GENES)
import matplotlib.patches as mpatches
from collections import Counter


import os
from collections import Counter

def detect_multiple_samples(raw_data_dir, list_name):
    """
    Detect sample IDs from filenames in a directory.
    Uses filename stems to infer sample IDs.
    """
    sample_names = []

    for fname in os.listdir(raw_data_dir):
        fpath = os.path.join(raw_data_dir, fname)
        
        
        if not os.path.isfile(fpath):
            continue

        name, ext = os.path.splitext(fname)
        if ext not in [".tsv", ".csv", ".mtx", ".h5", ".h5ad",".txt"]:  
            continue  # skip junk files

        # handle double extensions like ".counts.tsv"
        if ext in [".tsv", ".csv", ".mtx", ".h5",".txt"]:
            name = name.split(".")[0]

        # try to parse ID
        parts = name.split("_")
        if len(parts) > 1 and parts[1].startswith("T"):
            sample_id = parts[1]
        else:
            sample_id = name

        for sub in ["-barcodes", "-features", "-matrix", "_gene_expression"]:
            if sub in sample_id:
                sample_id = sample_id.split(sub)[0]

        if sample_id.startswith("GSM"):
            sample_id = "_".join(sample_id.split("_")[1:])

        

        sample_names.append(sample_id)

    ALL_SAMPLES = sorted(set(sample_names))

    print(f"[{raw_data_dir}] detected samples:", ALL_SAMPLES)

    add_list_to_config(
        list_name=list_name,
        list_values=ALL_SAMPLES,
        config_path="python_scripts/config.py"
    )

    return ALL_SAMPLES

def detect_replicates(raw_data_dir, base_list_name):
    """
    Detect base sample IDs from filenames.
    Example:
      F_14W_embryo1_1.h5ad → F_14W
      F_14W_embryo2_3.h5ad → F_14W
    """
    base_names = []

    for fname in os.listdir(raw_data_dir):
        if not fname.endswith(".h5ad"):
            continue

        name = os.path.splitext(fname)[0]
        parts = name.split("_")

        if len(parts) >= 2:
            base_id = "_".join(parts[:2])  # e.g. F_14W
        else:
            base_id = name

        base_names.append(base_id)

    ALL_BASES = sorted(set(base_names))
    print(f"[{raw_data_dir}] detected merged base samples:", ALL_BASES)

    add_list_to_config(
        list_name=base_list_name,
        list_values=ALL_BASES,
        config_path="python_scripts/config.py"
    )
    return ALL_BASES


def detect_compressed_files(raw_data_dir):
    for fname in os.listdir(raw_data_dir):
        full = os.path.join(raw_data_dir, fname)
        if fname.endswith('.gz'):
            with gzip.open(full, 'rb') as f_in:
                with open(full[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(full)  # Remove the original .gz file
###############################################################################################

import os
import scanpy as sc
import anndata as ad

def combine_replicates(base_sample, raw_dir):
    """
    Combine all .h5ad files that share the same base sample (e.g. F_14W).
    Preserves embryo ID and replicate ID in adata.obs.
    """

    files = [f for f in os.listdir(raw_dir) 
             if f.startswith(base_sample) and f.endswith(".h5ad")]

    if not files:
        raise FileNotFoundError(f"No replicates found for {base_sample} in {raw_dir}")

    adatas = []
    for f in files:
        a = sc.read(os.path.join(raw_dir, f))

        name = os.path.splitext(f)[0]
        parts = name.split("_")

        # defaults
        embryo, rep = None, None

        # find "embryoX" if present
        for p in parts:
            if p.startswith("embryo"):
                embryo = p.replace("embryo", "")
                break

        # take last token as replicate if it’s a digit
        if parts[-1].isdigit():
            rep = parts[-1]

        a.obs["sample"] = base_sample       # e.g. "F_14W"
        if embryo is not None:
            a.obs["embryo"] = embryo        # "1", "2", etc.
        if rep is not None:
            a.obs["replicate"] = rep        # "1", "2", etc.

        adatas.append(a)

    return ad.concat(adatas, join="outer")



def import_raw_data(sample, raw_file_type, raw_data_dir):
    """
    Expects files (for .mtx):
      *_matrix_{sample}.mtx(.gz)
      *_features_{sample}.tsv(.gz)
      *_barcodes_{sample}.tsv(.gz)
    Returns: adata_init, adata_init_path
    """
    adata_init = None
    adata_init_path = os.path.join(raw_data_dir, f"{sample}.h5ad")

    if raw_file_type == ".mtx":
        counts = obs = var = None
        for fname in os.listdir(raw_data_dir):
            if sample not in fname:
                continue
            full = os.path.join(raw_data_dir, fname)
            if "matrix" in fname:
                counts = csr_matrix(mmread(full)).T

            elif "features" in fname or "genes" in fname:
                var = pd.read_csv(full, header=None, sep="\t")

                if var.shape[1] >= 3:
                    var = var.iloc[:, :3]
                    var.columns = ["ensembl_id", "gene_id", "feature_type"]
                elif var.shape[1] == 2:
                    var.columns = ["ensembl_id", "gene_id"]
                else:
                    var.columns = ["ensembl_id"]
                    var["gene_id"] = var["ensembl_id"]
            elif "barcodes" in fname:
                obs = pd.read_csv(full, header=None, index_col=0)

                obs.index.name = None

        if counts is None or var is None or obs is None:
            raise FileNotFoundError(f"[{sample}] Missing one of counts/rows/cols files in {raw_data_dir}")

        adata_init = sc.AnnData(X=counts, obs=obs, var=var)
        
        
    elif raw_file_type == ".h5":
        target = None
        for fname in os.listdir(raw_data_dir):
            if sample in fname and fname.endswith(".h5"):
                target = os.path.join(raw_data_dir, fname)
                break
        if target is None:
            raise FileNotFoundError(f"[{sample}] .h5 file not found in {raw_data_dir}")

        # For 10x Genomics HDF5 files (CellRanger output), use sc.read_10x_h5
        adata_init = sc.read_10x_h5(target)

        # Make sure gene symbols are unique
        adata_init.var_names_make_unique()

    elif raw_file_type in (".csv", ".tsv", ".txt"):
        umi_data = None
        for fname in os.listdir(raw_data_dir):
            if sample in fname and ("counts" in fname or "gene_expression" in fname):
                full = os.path.join(raw_data_dir, fname)
                if fname.endswith(".csv") or ".csv" in fname:
                    umi_data = pd.read_csv(full, index_col=0)
                else:
                    umi_data = pd.read_csv(full, index_col=0, sep="\t")
                break
        if umi_data is None:
            raise FileNotFoundError(f"[{sample}] counts file not found in {raw_data_dir}")
        adata_init = sc.AnnData(X=umi_data)
        adata_init.var_names_make_unique()
        adata_init.raw = adata_init

    elif raw_file_type == ".h5ad":
        target = None
        for fname in os.listdir(raw_data_dir):
            if sample in fname and fname.endswith(".h5ad"):
                target = os.path.join(raw_data_dir, fname)
                break
        if target is None:
            raise FileNotFoundError(f"[{sample}] .h5ad file not found in {raw_data_dir}")
        adata_init = sc.read_h5ad(target)

    else:
        raise ValueError(f"Unsupported raw_file_type: {raw_file_type}")

    # --- now apply swap/cleanup to all cases ---
    import anndata

    print("First var name:", adata_init.var_names[0:10])
    print("First obs name:", adata_init.obs_names[0:10])

    first_var = str(adata_init.var_names[0])

    looks_like_barcode = (
        first_var.startswith(("AAA", "BT", "Fim"))
        or first_var.isupper() and first_var.startswith(("AAAC", "AAAG", "AAAT"))
        or len(first_var) > 12  # barcodes are usually long, genes usually short
)
    if looks_like_barcode:
        print("⚠️ obs/var appear swapped, rebuilding AnnData...")
        X = adata_init.X.T
        obs = adata_init.var.copy()
        var = adata_init.obs.copy()
        adata_fixed = anndata.AnnData(X=X, obs=obs, var=var)
        for layer in adata_init.layers:
            adata_fixed.layers[layer] = adata_init.layers[layer].T
        adata_init = adata_fixed
        print("fixed: cells =", adata_init.n_obs, "| genes =", adata_init.n_vars)
    else:
        print("No swap needed.")

    # --- prefer HGNC symbols ---
    if "feature_name" in adata_init.var.columns:        # NEW
        adata_init.var_names = adata_init.var["feature_name"].astype(str)  # NEW
    elif "gene_id" in adata_init.var.columns:
        adata_init.var_names = adata_init.var["gene_id"].astype(str)
    elif "ensembl_id" in adata_init.var.columns:
        adata_init.var_names = adata_init.var["ensembl_id"].astype(str)
    adata_init.var_names_make_unique()

    # --- resolve clashes ---
    if "gene_id" in adata_init.var.columns:
        if adata_init.var["gene_id"].astype(str).equals(
            pd.Index(adata_init.var_names).astype(str)
        ):
            adata_init.var.drop(columns=["gene_id"], inplace=True)
        else:
            adata_init.var.rename(columns={"gene_id": "hgnc_symbol"}, inplace=True)

        # --- tidy indices ---
    adata_init.var.index.name = None
    adata_init.obs.index.name = None

    # keep a raw copy
    adata_init.raw = adata_init.copy()
    adata_init.raw.var.index.name = None

    # --- save and return ---
    adata_init.write(adata_init_path)
    return adata_init, adata_init_path



############################################################################################

import scanpy as sc
import anndata as ad
import os

def combine_replicates(base_sample, raw_dir):
    # find all files like F_14W_embryo1_*.h5ad
    files = [f for f in os.listdir(raw_dir) if f.startswith(base_sample + "_") and f.endswith(".h5ad")]

    if not files:
        # fallback to single file
        return sc.read(os.path.join(raw_dir, f"{base_sample}.h5ad"))

    # load and tag replicates
    adatas = []
    for f in files:
        a = sc.read(os.path.join(raw_dir, f))
        rep_id = f.replace(base_sample + "_", "").replace(".h5ad", "")
        a.obs["sample"] = base_sample
        a.obs["replicate"] = rep_id
        adatas.append(a)

    return ad.concat(adatas, join="outer")


###############################################################################################

def prep_adata(adata, sample, figures_dir, tables_dir, 
    # Basic filtering
    min_genes=200,
    min_cells=10,
    # QC thresholds
    max_genes=2500,
    max_mt=5,
    # Outlier MAD thresholds
    n_mads_counts=5,
    n_mads_genes=5,
    n_mads_mt=3,
    n_mads_ribo=5,
    # Normalization
    normalize_target=1e4,
    save_log_norm=True,
    # HVGs
    n_top_hvgs=None,
    # Scaling
    scale_max=10,
    # Gene prefixes
    mt_prefix="MT-",
    ribo_prefix=("RPS", "RPL"),
):
    """
    Standard single-cell preprocessing/QC + UMAP on an AnnData object.
    """

    # ---- Raw counts → layers['counts']
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.raw.X.copy() if (adata.raw is not None and adata.raw.X is not None) else adata.X.copy()

    Xc = adata.layers["counts"]

    print("First 10 genes:", adata.var_names[:10].tolist())

    if adata.X.shape[1] != adata.var.shape[0]:
        print(f"⚠️ Mismatch in features: X has {adata.X.shape[1]} columns, var has {adata.var.shape[0]} rows.")
        adata.var = adata.var.iloc[:adata.X.shape[1], :].copy()
        adata.var_names = adata.var_names.astype(str)
    
    # --- Basic filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # --- Annotate MT + ribo
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(ribo_prefix)

    # --- QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True, percent_top=[20]
    )

    # --- Hard thresholds
    before = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :].copy()
    adata = adata[adata.obs.pct_counts_mt < max_mt, :].copy()
    after = adata.n_obs
    print(f"Hard QC filter: {before} → {after} cells retained")

    # --- MAD-based outlier detection
    def is_outlier(vec, nmads):
        if nmads is None:
            return np.zeros_like(vec, dtype=bool)  # no outliers if threshold disabled
        return (vec < np.median(vec) - nmads * median_abs_deviation(vec)) | (
            vec > np.median(vec) + nmads * median_abs_deviation(vec))


    adata.obs["outlier"] = (
        is_outlier(adata.obs["log1p_total_counts"], n_mads_counts)
        | is_outlier(adata.obs["log1p_n_genes_by_counts"], n_mads_genes)
        | is_outlier(adata.obs["pct_counts_mt"], n_mads_mt)
        | is_outlier(adata.obs["pct_counts_ribo"], n_mads_ribo)
    )

    before = adata.n_obs
    adata = adata[~adata.obs["outlier"], :].copy()
    after = adata.n_obs
    print(f"MAD filter: {before} → {after} cells retained")

    # ---- after QC filtering ----
    if adata.n_obs == 0:
        print(f"[{sample}] No cells left after QC. Skipping this sample.")
        return adata  # or return None if you prefer to drop it


    # --- Normalize + log transform
    sc.pp.normalize_total(adata, target_sum=normalize_target)
    sc.pp.log1p(adata)

    if save_log_norm:
        adata.layers["log_norm"] = adata.X.copy()

    # --- HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvgs)

    # Force-keep IMPORTANT_GENES + CTCFL/CTCF
    forced_genes = IMPORTANT_GENES + ["CTCFL", "CTCF"]
    forced_mask = adata.var_names.str.upper().isin([g.upper() for g in forced_genes])

    # Only keep forced genes that are actually expressed
    if sp.issparse(adata.X):
        nonzero_mask = np.array((adata.X > 0).sum(axis=0)).ravel() > 0
    else:
        nonzero_mask = (adata.X > 0).sum(axis=0) > 0

    forced_mask = forced_mask & nonzero_mask

    # Union HVGs + forced genes
    hvgs_or_forced = adata.var["highly_variable"] | forced_mask
    adata = adata[:, hvgs_or_forced].copy()

    print(f"Selected {adata.n_vars} genes (HVGs + forced genes).")
    missing = [g for g in forced_genes if g.upper() not in adata.var_names.str.upper()]
    if missing:
        print(f"⚠️ Skipped missing genes (not in dataset or all zero): {missing}")

    print(f"Selected {adata.n_vars} highly variable genes")

    # --- Regress out covariates
    #sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt", "pct_counts_ribo"])

    # --- Scale
    sc.pp.scale(adata, max_value=scale_max)

    # ---- PCA / neighbors / UMAP
    limit = min(adata.n_obs, adata.n_vars)
    if limit <= 1:
        print(f"[{sample}] Skipping PCA — not enough cells/genes (limit={limit})")
    else:
        n_comps = min(50, limit - 1)  # safe upper bound
        sc.pp.pca(adata, n_comps=n_comps, svd_solver="randomized")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_comps)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)




    # Get cell markers using helper function make_markers_df()
    
    if "celltype" in adata.obs.columns:
        adata.obs.rename(columns={"celltype": "cell_type"}, inplace=True)
       
    get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)

    # 9) CTCFL and CTCF expression - raw gene counts using helper function get_raw_gene_counts()
    print("Number of genes (processed):", adata.n_vars)
    print("First 10 genes (processed):", adata.var_names[:10].tolist())

    for gene in ["CTCF", "CTCFL"]:
        in_proc = gene in adata.var_names
        in_proc_case = any(g.upper() == gene for g in adata.var_names)
        in_counts = "counts" in adata.layers and gene in adata.var_names
        in_raw = adata.raw is not None and gene in adata.raw.var_names

        print(f"\n--- {gene} ---")
        print(f"Processed (var_names): {in_proc} (case-insensitive: {in_proc_case})")
        print(f"Counts layer: {in_counts}")
        print(f"Raw slot: {in_raw}")

        try:
            gene_expr = get_raw_gene_counts(adata, gene)
            total = gene_expr.sum()
            print(f"✅ {gene} total counts = {total} across {adata.n_obs} cells")
        except KeyError:
            print(f"❌ {gene} not found in any gene list")
            # fallback: mark all cells as negative
            adata.obs[f"{gene}_expr"] = f"{gene}_neg"
            adata.obs[f"{gene}_expr"] = adata.obs[f"{gene}_expr"].astype("category")
            continue  # skip export

        # classify cells
        pos_mask = gene_expr > 0
        adata.obs[f"{gene}_expr"] = np.where(pos_mask, f"{gene}_pos", f"{gene}_neg")
        adata.obs[f"{gene}_expr"] = adata.obs[f"{gene}_expr"].astype("category")

        print(f"{pos_mask.sum()} / {adata.n_obs} cells {gene}_pos")

        # only export CSVs for CTCFL
        if gene == "CTCFL":
            # choose raw-like source
            if adata.raw is not None:
                X_raw, var_names = adata.raw.X, adata.raw.var_names
            elif "counts" in adata.layers:
                X_raw, var_names = adata.layers["counts"], adata.var_names
            else:
                X_raw, var_names = adata.X, adata.var_names

            X_pos = X_raw[pos_mask, :]
            X_neg = X_raw[~pos_mask, :]

            if sp.issparse(X_pos):
                X_pos = X_pos.toarray()
            if sp.issparse(X_neg):
                X_neg = X_neg.toarray()

            df_pos = pd.DataFrame(X_pos, index=adata.obs_names[pos_mask], columns=var_names)
            df_neg = pd.DataFrame(X_neg, index=adata.obs_names[~pos_mask], columns=var_names)

            outdir = os.path.join(tables_dir, "gene")
            os.makedirs(outdir, exist_ok=True)

            df_pos.to_csv(os.path.join(outdir, f"{sample}_{gene}_pos.csv"))
            #df_neg.to_csv(os.path.join(outdir, f"{sample}_{gene}_neg.csv"))

            print(f"{sample}: {df_pos.shape[0]} {gene}_pos cells and {df_neg.shape[0]} {gene}_neg cells")




    # ---- Save UMAP
    global CELL_TYPE_COLORS
    global GENE_COLOR_MAP
    if "CELL_TYPE_COLORS" not in globals():
        CELL_TYPE_COLORS = {}  # fallback empty map

    def safe_palette(adata, column, base_palette):
        cats = adata.obs[column].cat.categories if column in adata.obs else []
        return {c: base_palette.get(c, "#d3d3d3") for c in cats}
    

    os.makedirs(os.path.join(figures_dir, "umap"), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(18,6))

    # Panel A: cell types
    sc.pl.umap(
        adata,
        color="cell_type",
        title=f"{sample} cell types",
        frameon=False,
        legend_fontweight="normal",
        legend_fontsize=10,
        show=False,
        ax=axes[0],
        palette=CELL_TYPE_COLORS,   # <-- use universal cell type map
    )

    # Panel B: CTCF status
    sc.pl.umap(
        adata,
        color="CTCF_expr",
        title=f"{sample} CTCF status",
        frameon=False,
        legend_fontweight="normal",
        legend_fontsize=10,
        show=False,
        ax=axes[1],
        palette=safe_palette(adata, "CTCF_expr", GENE_COLOR_MAP),
    )

    # Panel C: CTCFL status
    sc.pl.umap(
        adata,
        color="CTCFL_expr",
        title=f"{sample} CTCFL status",
        frameon=False,
        legend_fontweight="normal",
        legend_fontsize=10,
        show=False,
        ax=axes[2],
        palette=safe_palette(adata, "CTCFL_expr", GENE_COLOR_MAP),
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(figures_dir, "umap", f"{sample}_umap_panels.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


    print("CTCFL kept?", "CTCFL" in adata.var_names)
    return adata


def make_adata(
    adata_init_path,
    adata_working_path,
    sample,
    figures_dir,
    tables_dir,
    # --- QC params with defaults ---
    min_genes=200,
    min_cells=10,
    max_genes=2500,
    max_mt=5,
    n_mads_counts=5,
    n_mads_genes=5,
    n_mads_mt=3,
    n_mads_ribo=5,
    normalize_target=1e4,
    n_top_hvgs=None,
    scale_max=10,
    mt_prefix="MT-",
    ribo_prefix=("RPS", "RPL"),
):
    print(adata_init_path, adata_working_path)

    if os.path.exists(adata_working_path):
        adata = sc.read_h5ad(adata_working_path)
    else:
        adata = sc.read_h5ad(adata_init_path)

        # ✅ pass through the params, don’t hardcode them
        adata = prep_adata(
            adata,
            sample,
            figures_dir,
            tables_dir,
            min_genes=min_genes,
            min_cells=min_cells,
            max_genes=max_genes,
            max_mt=max_mt,
            n_mads_counts=n_mads_counts,
            n_mads_genes=n_mads_genes,
            n_mads_mt=n_mads_mt,
            n_mads_ribo=n_mads_ribo,
            normalize_target=normalize_target,
            n_top_hvgs=n_top_hvgs,
            scale_max=scale_max,
            mt_prefix=mt_prefix,
            ribo_prefix=ribo_prefix,
        )

        print(f"Saving processed data to {adata_working_path}")
        adata.write(adata_working_path)

    return adata


################################################################################################
def make_bdata(donor, adata_init_path, working_adata_path, working_bdata_file_path, figures_dir, tables_dir):
    """
    Subset a processed AnnData to a specific donor/sample and persist the subset.

    Parameters
    ----------
    adata : sc.AnnData
        Preprocessed AnnData containing multiple donors/samples.
    sample : str
        Value in adata.obs['donor_id'] to subset.
    working_bdata_file_path : str
        Output path for the subset AnnData (.h5ad).

    Returns
    -------
    sc.AnnData

    Usage
    -----
    bdata = make_bdata(adata, sample, working_bdata_file_path)
    """
    if os.path.exists(working_bdata_file_path):
        # Fast path: read the subset if it already exists
        bdata = sc.read_h5ad(working_bdata_file_path)

    else:

        if os.path.exists(working_adata_path):
            # Reuse previously processed object
            adata = sc.read_h5ad(working_adata_path)
        else:
            # Read raw and preprocess, then persist
            adata = sc.read_h5ad(adata_init_path)

   
        # Subset by donor_id and persist
        bdata = adata[adata.obs["donor_id"] == donor].copy()

        bdata = prep_adata(bdata, donor, figures_dir, tables_dir)
        print(f"Saving processed data to {working_bdata_file_path}")
        #os.makedirs(os.path.dirname(adata_working_path), exist_ok=True)
        bdata.write(working_bdata_file_path)
        os.makedirs(os.path.dirname(working_bdata_file_path), exist_ok=True)
        print(f"Saving processed data to {working_bdata_file_path}")
        bdata.write(working_bdata_file_path)

    return bdata
