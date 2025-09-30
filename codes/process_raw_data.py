# codes/process_adata.py

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

#### Custom ####

from setup import add_list_to_config, GENE_COLOR_MAP, CELL_TYPE_COLORS


#### Functions ####

def detect_multiple_samples(raw_data_dir, list_name):
    """
    This function is for detecting sample IDs (not replicates) from filenames in a directory.
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
        config_path="config.py"
    )

    return ALL_SAMPLES

def detect_replicates(raw_data_dir, base_list_name):
    """
    Detect replicates (not base samples) only after processing raw data files to .h5ad.
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
        config_path="config.py"
    )
    return ALL_BASES


def import_raw_data(sample, raw_data_dir, raw_file_type):
    """
    Import raw single-cell data for a given sample.

    Returns:
        adata_init (AnnData), adata_init_path (str)
    """
    adata_init = None
    adata_init_path = os.path.join(raw_data_dir, f"{sample}.h5ad")

    
    
    for fname in os.listdir(raw_data_dir):
        if sample not in fname:
            continue

   

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
        target = next((os.path.join(raw_data_dir, f) for f in os.listdir(raw_data_dir)
                       if sample in f and f.endswith(".h5")), None)
        if target is None:
            raise FileNotFoundError(f"[{sample}] .h5 file not found in {raw_data_dir}")
        adata_init = sc.read_10x_h5(target)
        adata_init.var_names_make_unique()

    elif raw_file_type in (".csv", ".tsv", ".txt"):
        umi_data = None
        for fname in os.listdir(raw_data_dir):
            if sample in fname and ("counts" in fname or "gene_expression" in fname):
                full = os.path.join(raw_data_dir, fname)
                sep = "," if fname.endswith(".csv") or ".csv" in fname else "\t"
                umi_data = pd.read_csv(full, index_col=0, sep=sep)
                break
        if umi_data is None:
            raise FileNotFoundError(f"[{sample}] counts file not found in {raw_data_dir}")
        adata_init = sc.AnnData(X=umi_data)
        adata_init.var_names_make_unique()
        adata_init.raw = adata_init

    elif raw_file_type == ".h5ad":
        target = next((os.path.join(raw_data_dir, f) for f in os.listdir(raw_data_dir)
                       if sample in f and f.endswith(".h5ad")), None)
        if target is None:
            raise FileNotFoundError(f"[{sample}] .h5ad file not found in {raw_data_dir}")
        adata_init = sc.read_h5ad(target)

    else:
        raise ValueError(f"Unsupported raw_file_type: {raw_file_type}")

    # --- check for swapped obs/var ---
    print("First var name:", adata_init.var_names[:10])
    print("First obs name:", adata_init.obs_names[:10])
    first_var = str(adata_init.var_names[0])

    looks_like_barcode = (
        first_var.startswith(("AAA", "BT", "Fim"))
        or (first_var.isupper() and first_var.startswith(("AAAC", "AAAG", "AAAT")))
        or len(first_var) > 12
    )
    if looks_like_barcode:
        print("[WARN] obs/var appear swapped, rebuilding AnnData...")
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

    # --- set var_names preference ---
    if "feature_name" in adata_init.var.columns:
        adata_init.var_names = adata_init.var["feature_name"].astype(str)
    elif "gene_id" in adata_init.var.columns:
        adata_init.var_names = adata_init.var["gene_id"].astype(str)
    elif "ensembl_id" in adata_init.var.columns:
        adata_init.var_names = adata_init.var["ensembl_id"].astype(str)
    adata_init.var_names_make_unique()

    # --- resolve clashes ---
    if "gene_id" in adata_init.var.columns:
        if adata_init.var["gene_id"].astype(str).equals(pd.Index(adata_init.var_names).astype(str)):
            adata_init.var.drop(columns=["gene_id"], inplace=True)
        else:
            adata_init.var.rename(columns={"gene_id": "hgnc_symbol"}, inplace=True)

    # --- tidy indices ---
    adata_init.var.index.name = None
    adata_init.obs.index.name = None

    # keep a raw copy
    adata_init.raw = adata_init.copy()
    adata_init.raw.var.index.name = None

    # save and return
    adata_init.write(adata_init_path)
    return adata_init, adata_init_path

######## Process adata #########

def prep_adata(
    adata,
    sample,
    figures_dir,
    tables_dir,
    # --- QC params ---
    min_genes=200,
    min_cells=10,
    max_genes=2500,
    max_mt=5,
    n_mads_counts=5,
    n_mads_genes=5,
    n_mads_mt=3,
    n_mads_ribo=5,
    # --- Normalization ---
    normalize_target=1e4,
    save_log_norm=True,
    # --- HVGs ---
    n_top_hvgs=None,
    forced_genes=None,          # NEW: pass your own list, instead of relying on globals
    # --- Scaling ---
    scale_max=10,
    # --- Gene prefixes ---
    mt_prefix="MT-",
    ribo_prefix=("RPS", "RPL"),
    # --- Color maps ---
    gene_color_map=None,
    cell_type_colors=None,
    # --- External helpers ---
    markers_file_path=None,
    get_cell_type_fn=None,
    get_raw_gene_counts_fn=None,
):
    """
    Standard single-cell preprocessing/QC + UMAP on an AnnData object.
    """

    # ---- Raw counts → layers['counts']
    if "counts" not in adata.layers:
        adata.layers["counts"] = (
            adata.raw.X.copy() if (adata.raw is not None and adata.raw.X is not None) else adata.X.copy()
        )

    print("First 10 genes:", adata.var_names[:10].tolist())

    # ---- Safety check for var alignment
    if adata.X.shape[1] != adata.var.shape[0]:
        print(f"[WARN] Mismatch: X has {adata.X.shape[1]} cols, var has {adata.var.shape[0]} rows")
        adata.var = adata.var.iloc[:adata.X.shape[1], :].copy()
        adata.var_names = adata.var_names.astype(str)

    # ---- Basic filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # ---- Annotate MT + ribo
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(ribo_prefix)

    # ---- QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True, percent_top=[20])

    # ---- Hard thresholds
    before = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :].copy()
    adata = adata[adata.obs.pct_counts_mt < max_mt, :].copy()
    after = adata.n_obs
    print(f"Hard QC filter: {before} → {after} cells retained")

    # ---- MAD-based outlier detection
    def is_outlier(vec, nmads):
        if nmads is None:
            return np.zeros_like(vec, dtype=bool)
        med = np.median(vec)
        mad = median_abs_deviation(vec)
        return (vec < med - nmads * mad) | (vec > med + nmads * mad)

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

    if adata.n_obs == 0:
        print(f"[{sample}] No cells left after QC. Skipping.")
        return adata

    # ---- Normalize + log
    sc.pp.normalize_total(adata, target_sum=normalize_target)
    sc.pp.log1p(adata)
    if save_log_norm:
        adata.layers["log_norm"] = adata.X.copy()

    # ---- HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvgs)

    # ---- Forced genes
    if forced_genes is None:
        forced_genes = ["CTCFL", "CTCF"]  # fallback default
    forced_mask = adata.var_names.str.upper().isin([g.upper() for g in forced_genes])

    # keep only expressed forced genes
    if sp.issparse(adata.X):
        nonzero_mask = np.array((adata.X > 0).sum(axis=0)).ravel() > 0
    else:
        nonzero_mask = (adata.X > 0).sum(axis=0) > 0
    forced_mask = forced_mask & nonzero_mask

    hvgs_or_forced = adata.var["highly_variable"] | forced_mask
    adata = adata[:, hvgs_or_forced].copy()
    print(f"Selected {adata.n_vars} genes (HVGs + forced).")

    missing = [g for g in forced_genes if g.upper() not in adata.var_names.str.upper()]
    if missing:
        print(f"[WARN] Missing forced genes: {missing}")

    # ---- Scale
    sc.pp.scale(adata, max_value=scale_max)

    # ---- PCA / neighbors / UMAP
    limit = min(adata.n_obs, adata.n_vars)
    if limit > 1:
        n_comps = min(50, limit - 1)
        sc.pp.pca(adata, n_comps=n_comps, svd_solver="randomized")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_comps)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)

    # ---- Cell type inference (optional)
    if "celltype" in adata.obs.columns:
        adata.obs.rename(columns={"celltype": "cell_type"}, inplace=True)
    if get_cell_type_fn is not None and markers_file_path:
        get_cell_type_fn(adata, markers_file_path=markers_file_path)

    # ---- Gene-specific checks
    if get_raw_gene_counts_fn is not None:
        for gene in ["CTCF", "CTCFL"]:
            try:
                gene_expr = get_raw_gene_counts_fn(adata, gene)
                pos_mask = gene_expr > 0
                adata.obs[f"{gene}_expr"] = np.where(pos_mask, f"{gene}_pos", f"{gene}_neg")
                adata.obs[f"{gene}_expr"] = adata.obs[f"{gene}_expr"].astype("category")
                print(f"{sample}: {pos_mask.sum()} {gene}_pos cells")
            except Exception:
                print(f"{gene} not found. Marking all as negative.")
                adata.obs[f"{gene}_expr"] = f"{gene}_neg"
                adata.obs[f"{gene}_expr"] = adata.obs[f"{gene}_expr"].astype("category")

    # ---- Plot UMAPs
    def safe_palette(column, base_palette):
        cats = adata.obs[column].cat.categories if column in adata.obs else []
        return {c: base_palette.get(c, "#d3d3d3") for c in cats}

    os.makedirs(os.path.join(figures_dir, "umap"), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))

    sc.pl.umap(
    adata, color="cell_type", title=f"{sample} cell types",
    frameon=False, legend_fontsize=10, show=False, ax=axes[0],
    palette=cell_type_colors if cell_type_colors is not None else CELL_TYPE_COLORS
    )

    sc.pl.umap(
        adata, color="CTCFL_expr", title=f"{sample} CTCFL",
        frameon=False, legend_fontsize=10, show=False, ax=axes[1],
        palette=safe_palette("CTCFL_expr", gene_color_map if gene_color_map is not None else GENE_COLOR_MAP)
    )

    sc.pl.umap(
        adata, color="CTCF_expr", title=f"{sample} CTCF",
        frameon=False, legend_fontsize=10, show=False, ax=axes[2],
        palette=safe_palette("CTCF_expr", gene_color_map if gene_color_map is not None else GENE_COLOR_MAP)
    )


    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "umap", f"{sample}_umap_panels.png"), dpi=300)
    plt.close()

    return adata

##### make adata #####

def make_adata(
    adata_init_path,
    adata_working_path,
    sample,
    figures_dir,
    tables_dir,
    # --- QC params ---
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
    # --- Optional customizations ---
    forced_genes=None,
    gene_color_map=None,
    cell_type_colors=None,
    markers_file_path=None,
    get_cell_type_fn=None,
    get_raw_gene_counts_fn=None,
):
    """
    Wrapper around prep_adata to build a processed AnnData.
    Reuses adata_working_path if it exists, otherwise processes fresh.
    """
    print(adata_init_path, adata_working_path)

    if os.path.exists(adata_working_path):
        adata = sc.read_h5ad(adata_working_path)
    else:
        adata = sc.read_h5ad(adata_init_path)

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
            forced_genes=forced_genes,
            gene_color_map=gene_color_map,
            cell_type_colors=cell_type_colors,
            markers_file_path=markers_file_path,
            get_cell_type_fn=get_cell_type_fn,
            get_raw_gene_counts_fn=get_raw_gene_counts_fn,
        )

        print(f"Saving processed data to {adata_working_path}")
        adata.write(adata_working_path)

    return adata
