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

import decoupler as dc

# Plot defaults
sc.set_figure_params(figsize=(3, 3), frameon=False)

# Project configuration (paths, constants, etc.)
from config import *


def add_list_to_config(list_name, list_values, config_path="python_scripts/config.py"):

    # read file
    with open(config_path, "r") as f:
        lines = f.readlines()

    found = False
    with open(config_path, "w") as f:
        for line in lines:
            if line.strip().startswith(f"{list_name} ="):
                f.write(f"{list_name} = {list_values}\n")
                found = True
            else:
                f.write(line)
        if not found:
            f.write(f"\n{list_name} = {list_values}\n")
    print(f"set {list_name} in config.py")


def gene_name_to_ensembl(adata):
    #Make the gene_id as the gene index
    adata.var["hgnc_symbol"] = adata.var["gene_id"].astype(str)
    adata.var = adata.var[["ensembl_id", "hgnc_symbol"]] #preserve onlyy hgnc and ensembl id
    adata.var_names = adata.var["hgnc_symbol"] #Use symbols as var_names
    adata.var_names_make_unique()
    adata.var.index.name = None #make sure the index name does NOT clash with a column name

    #Freeze counts into a permanent raw layer
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy() #stores frozen snapshot of counts and metadata for plotting marker genes after normalizaion

    # Check
    print(adata)
    print("var head:\n", adata.var.head())

        # ---- Optional: adopt 'feature_name' as var_names if available
    if "feature_name" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names.astype(str)
        adata.var_names = adata.var["feature_name"].astype(str)
        # avoid write-time collisions
        adata.var.drop(columns=["feature_name"], inplace=True)
    adata.var.index.name = None
    if adata.raw is not None and adata.raw.var is not None:
        adata.raw.var.index.name = None


    return adata


def make_markers_df(markers_file_path = MARKERS_FILE_PATH):
    """Usage: markers_df = make_markers_df()"""
    markers_file_path = "/vf/users/scottaa/human_fetal_gonad_rnaseq/results/tables/markers_df.csv"
    if os.path.exists(markers_file_path):
        markers_df = pd.read_csv(markers_file_path, index_col = 0)
    else:

        markers = dc.op.resource("PanglaoDB", organism="human")
        # Filter by canonical_marker and human
        markers = markers[
        markers["human"].astype(bool)
        & markers["canonical_marker"].astype(bool)
        & (markers["human_sensitivity"].astype(float) > 0.5)
        ]

        # Remove duplicated entries
        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

        # Format
        markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
        markers = markers[["source", "target"]]

        markers_df = pd.DataFrame(markers)

        print(markers_df)

        markers_df.to_csv(os.path.join("/vf/users/scottaa/human_fetal_gonad_rnaseq/results/tables","markers_df.csv"))
    
    return markers_df





def get_raw_gene_counts(adata, gene):
    if adata.raw is None:
        raise ValueError("adata.raw is None — no raw counts stored.")

    if gene not in adata.raw.var_names:
        raise KeyError(f"{gene} not found in adata.raw.var_names")

    # extract column directly
    Xg = adata.raw[:, gene].X
    return Xg.toarray().ravel() if sp.issparse(Xg) else np.ravel(Xg)
def get_cell_type(adata, markers_file_path=MARKERS_FILE_PATH):
    """
    Compute cell-type enrichment (ULM) and annotate:
      - adata.obs['predicted_cell_type'] (per-cell max score)
      - adata.obs['cell_type'] (Leiden cluster → top cell type)
    """

    # ---- load or fetch marker network (source-target long df) ----
    if os.path.exists(markers_file_path):
        markers_df = pd.read_csv(markers_file_path, index_col=0)
        # ensure it's only the two required cols
        if not set(["source", "target"]).issubset(markers_df.columns):
            raise ValueError("markers_file must have columns: ['source','target']")
        markers_df = markers_df[["source", "target"]]
    else:
        # decoupler ≥2.x
        try:
            markers = dc.get_resource(name="PanglaoDB", organism="human", license="academic")
        except AttributeError:
            # fallback for older builds
            markers = dc.op.resource("PanglaoDB", organism="human", license="academic")

        markers = markers[
            markers["human"].astype(bool)
            & markers["canonical_marker"].astype(bool)
            & (markers["human_sensitivity"].astype(float) > 0.5)
        ].copy()

        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
        markers_df = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})[["source", "target"]]
        markers_df.to_csv(markers_file_path)
    
    # ---- run ULM and get per-cell best label ----
    # ---- run decoupler with any available method, then read scores ----
    dc.mt.ulm(data=adata, net=markers_df, tmin=3)
    adata.obsm["score_ulm"]
    score = dc.pp.get_obsm(adata, key="score_ulm")
    
    # convert to DataFrame if it's an AnnData
    score_df = score.to_df() if hasattr(score, "to_df") else score

    # per-cell best label
    adata.obs["predicted_cell_type"] = score_df.idxmax(axis=1)

    score = dc.pp.get_obsm(adata, key="score_ulm")  # AnnData or DataFrame depending on version
    score_df = score.to_df() if hasattr(score, "to_df") else score
    adata.obs["predicted_cell_type"] = score_df.idxmax(axis=1)

    # ---- rank top cell types per Leiden cluster ----
    df = dc.tl.rankby_group(adata=score, groupby="leiden", reference="rest", method="t-test_overestim_var")
    df = df[df["stat"] > 0]

    n_ctypes = 3
    ctypes_dict = (df.groupby("group").head(n_ctypes)
                     .groupby("group")["name"].apply(list).to_dict())

    dict_ann = (df.groupby("group").head(1)
                  .set_index("group")["name"].to_dict())

    s = adata.obs["leiden"].astype(str)
    adata.obs["cell_type"] = pd.Categorical(s.map(dict_ann).fillna(s))

    return adata

def subset_adata(adata, sample, tables_dir=TABLES_DIR):
    outdir = os.path.join(tables_dir, "single_cell_analysis", sample)
    os.makedirs(outdir, exist_ok=True)

    # 1) Threshold on RAW only (so you get 1–2 positives again)
    ctcfl_vals = get_raw_gene_counts(adata, "CTCFL")
    pos_mask = ctcfl_vals > 0
    neg_mask = ~pos_mask

    # 2) Export from the SAME RAW layer
    if adata.raw is None:
        raise ValueError("adata.raw is None — cannot export raw count matrices.")
    X_raw = adata.raw.X
    var_names = adata.raw.var_names

    # slice first, THEN densify (memory-safe, same output as before)
    X_pos = X_raw[pos_mask, :]
    X_neg = X_raw[neg_mask, :]

    if sp.issparse(X_pos):
        X_pos = X_pos.toarray()
    else:
        X_pos = np.asarray(X_pos)

    if sp.issparse(X_neg):
        X_neg = X_neg.toarray()
    else:
        X_neg = np.asarray(X_neg)

    df_pos = pd.DataFrame(X_pos, index=adata.obs_names[pos_mask], columns=var_names)
    df_neg = pd.DataFrame(X_neg, index=adata.obs_names[neg_mask], columns=var_names)

    #df_pos.to_csv(os.path.join(outdir, f"{sample}_ctcfl_pos.csv"))
    #df_neg.to_csv(os.path.join(outdir, f"{sample}_ctcfl_neg.csv"))

    print(f"{sample}: {df_pos.shape[0]} CTCFL_pos cells and {df_neg.shape[0]} CTCFL_neg cells")
