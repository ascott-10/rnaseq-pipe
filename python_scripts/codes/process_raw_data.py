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

# Plot defaults
sc.set_figure_params(figsize=(3, 3), frameon=False)

# Project configuration (paths, constants, etc.)
from config import *
import gzip
import shutil
from codes.utils import make_markers_df, get_raw_gene_counts, get_cell_type, add_list_to_config

def detect_multiple_samples(raw_data_dir = RAW_DATA_DIR):
    sample_names = []
    mult_samples = []
    single_samples = []

    
    for fname in os.listdir(raw_data_dir):
        full = os.path.join(raw_data_dir, fname)
        sample_name = fname.split('_')
        if len(sample_name) > 1:
            sample_names.append(sample_name[1])

    for sample in sample_names:
        if sample_names.count(sample) > 1:
            mult_samples.append(sample)
            
        else:
            single_samples.append(sample)
            

    print("single samples:", single_samples)
    print("Multiple samples:", mult_samples)

    single_samples.sort()
    mult_samples.sort()

    ALL_SAMPLES = list(set(mult_samples))

    add_list_to_config(list_name = "ALL_SAMPLES", list_values = ALL_SAMPLES, config_path="python_scripts/config.py")
    
    return ALL_SAMPLES

def detect_compressed_files(raw_data_dir = RAW_DATA_DIR):
    for fname in os.listdir(raw_data_dir):
        full = os.path.join(raw_data_dir, fname)
        if fname.endswith('.gz'):
            with gzip.open(full, 'rb') as f_in:
                with open(full[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(full)  # Remove the original .gz file
###############################################################################################

def import_raw_data(sample, raw_file_type=RAW_FILE_TYPE, raw_data_dir=RAW_DATA_DIR):
    """
    Expects files (for .mtx):
      *_matrix_{sample}.mtx(.gz)
      *_features_{sample}.tsv(.gz)
      *_barcodes_{sample}.tsv(.gz)
    Returns: adata_init, adata_init_path
    """
    adata_init = None
    adata_init_path = os.path.join(RAW_DATA_DIR, f"{sample}.h5ad")

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
        # prefer HGNC symbols for var_names, fallback to Ensembl
        if "gene_id" in adata_init.var.columns:
            adata_init.var_names = adata_init.var["gene_id"].astype(str)
        elif "ensembl_id" in adata_init.var.columns:
            adata_init.var_names = adata_init.var["ensembl_id"].astype(str)
        adata_init.var_names_make_unique()

        # --- fix write clash: index name vs column name ---
        # If the column 'gene_id' is now the same as var_names, drop it; else keep it but rename.
        if "gene_id" in adata_init.var.columns:
            if adata_init.var["gene_id"].astype(str).equals(pd.Index(adata_init.var_names).astype(str)):
                adata_init.var.drop(columns=["gene_id"], inplace=True)
            else:
                adata_init.var.rename(columns={"gene_id": "hgnc_symbol"}, inplace=True)

        # Clear index name(s) to avoid conflicts when writing (also affects .raw)
        adata_init.var.index.name = None
        if adata_init.obs.index.name is not None:
            adata_init.obs.index.name = None

        # Set raw AFTER cleaning, then also clear raw.var index name
        adata_init.raw = adata_init.copy()
        adata_init.raw.var.index.name = None

    elif raw_file_type in (".csv", ".tsv"):
        umi_data = None
        for fname in os.listdir(raw_data_dir):
            if sample in fname and "counts" in fname:
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

    adata_init.write(adata_init_path)
    return adata_init, adata_init_path



###############################################################################################

def prep_adata(adata, sample, figures_dir = FIGURES_DIR):
    """
    Standard single-cell preprocessing/QC + UMAP on an AnnData object.

    Steps:
      - Store raw counts in a 'counts' layer (from adata.raw.X or adata.X)
      - Basic QC filtering (cells/genes; mito percent)
      - Normalize/log1p and store 'log_norm' layer
      - HVGs, regress out covariates, scale
      - PCA, neighbors, UMAP, Leiden
      - Save UMAP figure (figures_dir/umap/initial_umap.png)
      

    Returns
    -------
    sc.AnnData  (modified in place, also returned)
  """

    # ---- Raw counts → layers['counts']
    if "counts" not in adata.layers:
        if adata.raw is not None and adata.raw.X is not None:
            adata.layers["counts"] = adata.raw.X.copy()
        else:
            # create a raw snapshot if missing
            adata.raw = adata.copy()
            adata.layers["counts"] = adata.X.copy()

    # quick preview so you can sanity-check integers in logs
    Xc = adata.layers["counts"]
    _preview = Xc[:10, :10].toarray() if sp.issparse(Xc) else np.asarray(Xc[:10, :10])
    print(_preview)


    
    # ---- Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # ---- Mito + QC metrics
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # thresholds: tweak as needed
    adata = adata[adata.obs.n_genes_by_counts < 2500, :].copy()
    adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

    # ---- Normalize + log1p; keep a copy
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log_norm"] = adata.X.copy()

    # ---- HVGs / regress / scale
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    # ---- PCA / neighbors / UMAP
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)  # writes to adata.obsm['X_umap']

    # ---- Leiden clustering
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)


    # Get cell markers using helper function make_markers_df()
    if "cell_type" not in adata.obs.columns:
        if "celltype" in adata.obs.columns:
            adata.obs.rename(columns={"celltype": "cell_type"}, inplace=True)
        else:
            get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)

    #9)CTCFL and CTCF expression - raw gene counts using helper function get_raw_gene_counts()
    ctcf_expr  = get_raw_gene_counts(adata, 'CTCF')
    mask_ctcf = ctcf_expr > 0

    ctcfl_expr = get_raw_gene_counts(adata, 'CTCFL')
    mask_ctcfl = ctcfl_expr > 0
    
    # assign categories
    adata.obs["ctcf_expr"] = np.where(mask_ctcf, "CTCF_pos", "CTCF_neg")
    adata.obs["ctcfl_expr"] = np.where(mask_ctcfl, "CTCFL_pos", "CTCFL_neg")

    # make categorical
    adata.obs["ctcf_expr"] = adata.obs["ctcf_expr"].astype("category")
    adata.obs["ctcfl_expr"] = adata.obs["ctcfl_expr"].astype("category")

    print(f' \n {mask_ctcf.sum()} / {adata.n_obs} cells CTCF_pos')
    print(f'{mask_ctcfl.sum()} / {adata.n_obs} cells CTCFL_pos')
   
    print(f'adata.obs: {adata.obs}')


    # ---- Save UMAP
    os.makedirs(os.path.join(figures_dir, "umap"), exist_ok=True)
    sc.pl.umap(
        adata,
        color=["cell_type", "ctcfl_expr"],
        title=f'{sample} UMAP',
        frameon=False,
        legend_fontweight="normal",
        legend_fontsize=15,
        show=False,
    )
    plt.savefig(
        os.path.join(figures_dir, "umap", f'{sample}_initial_umap.png'),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    return adata


def make_adata(adata_init_path, adata_working_path, sample, figures_dir = FIGURES_DIR):
    """
    Load a preprocessed AnnData if available; otherwise read raw h5ad,
    run preprocessing (prep_adata), and save a processed copy.

    Parameters
    ----------
    adata_init_path : str
        Path to source AnnData (.h5ad).
    adata_working_path : str
        Path to processed AnnData (.h5ad) to read/write.
    figures_dir : str
        Directory to save figures generated during preprocessing.

    Returns
    -------
    sc.AnnData

    Usage
    -----
    adata = make_adata(adata_init_path, adata_working_path, figures_dir)
    """
    print(adata_init_path, adata_working_path)

    if os.path.exists(adata_working_path):
        # Reuse previously processed object
        adata = sc.read_h5ad(adata_working_path)
        
    else:
        # Read raw and preprocess, then persist
        adata = sc.read_h5ad(adata_init_path)
    
        adata = prep_adata(adata, sample, figures_dir = FIGURES_DIR)
        print(f"Saving processed data to {adata_working_path}")
        #os.makedirs(os.path.dirname(adata_working_path), exist_ok=True)
        adata.write(adata_working_path)

    return adata

################################################################################################
def make_bdata(adata, sample, working_bdata_file_path):
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
        # Subset by donor_id and persist
        bdata = adata[adata.obs["donor_id"] == sample].copy()
        os.makedirs(os.path.dirname(working_bdata_file_path), exist_ok=True)
        print(f"Saving processed data to {working_bdata_file_path}")
        bdata.write(working_bdata_file_path)

    return bdata
