# projects/E-MTAB-8559/python_scripts/codes/cell_annotation.py

import os
import numpy as np
import pandas as pd
import random

import scanpy as sc
import scipy.sparse as sp

import matplotlib.pyplot as plt
import decoupler as dc
import colorcet as cc

from config import *


# ---- globals ----
GENE_COLOR_MAP = {"None": "#d3d3d3"}
CELL_TYPE_COLORS = {}
GENE_COLOR_INDEX = 0

#### Make Gene List ####
def make_markers_df(markers_file_path=MARKERS_FILE_PATH):
    """
    Load or build a markers dataframe.
    Returns df with columns: source (cell type), target (gene).
    markers_df = make_markers_df(markers_file_path=MARKERS_FILE_PATH)
    """
    if os.path.exists(markers_file_path):
        markers_df = pd.read_csv(markers_file_path, index_col=0)
    else:
        markers = dc.op.resource("PanglaoDB", organism="mouse")
        markers = markers[
            markers["human"].astype(bool)
            & markers["canonical_marker"].astype(bool)
            & (markers["human_sensitivity"].astype(float) > 0.5)
        ]
        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
        markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
        markers = markers[["source", "target"]]

        markers_df = pd.DataFrame(markers)
        markers_df.to_csv(markers_file_path)

    return markers_df

def make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH):

    """CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)"""
    cell_markers_file_path = pd.read_csv(MARKERS_FILE_PATH)
    MARKER_GENES = list(cell_markers_file_path.target)
    dataset_file_path = pd.read_csv(ALL_GENES_FILE_PATH)
    dataset_cleaned_genes = []
    if "hgnc_symbol" in dataset_file_path.columns:
        dataset_genes =  list(dataset_file_path["hgnc_symbol"])
        dataset_cleaned_genes = pd.Series(dataset_genes).dropna().tolist()

    all_genes = list(set(MARKER_GENES + dataset_cleaned_genes))

    CLEANED_GENES = pd.Series(all_genes).dropna().tolist()

    return CLEANED_GENES

def make_color_map(items, fallback="#d3d3d3", palette="glasbey_light"):
    """
    Assign shuffled colors from a palette to items.
    Includes defaults for key BORIS genes and adds *_pos/*_neg variants.
    Global state: GENE_COLOR_MAP, GENE_COLOR_INDEX
    """
    global GENE_COLOR_INDEX

    # clean input: drop NaN, enforce string
    items = [str(g) for g in items if pd.notna(g) and str(g).strip() != ""]

    if palette in cc.palette:
        all_colors = list(cc.palette[palette])
    else:
        all_colors = list(getattr(cc, palette))

    random.seed(555)  # reproducible
    random.shuffle(all_colors)

    # ---- defaults first ----
    defaults = {
        "CTCFL":   "#04531c",
        "CTCF":    "#4b11d3",
    }
    GENE_COLOR_MAP.update(defaults)

    # ---- assign new items ----
    for g in items:
        if g not in GENE_COLOR_MAP:
            GENE_COLOR_MAP[g] = all_colors[GENE_COLOR_INDEX % len(all_colors)]
            GENE_COLOR_INDEX += 1

    # ---- enforce None ----
    GENE_COLOR_MAP["None"] = fallback

    # always return safe mapping
    return {g: GENE_COLOR_MAP.get(g, fallback) for g in items}


def init_color_maps(genes, markers_file_path=MARKERS_FILE_PATH, adata=None):
    """
    Initialize gene and cell-type color maps.

    GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(genes, markers_file_path=MARKERS_FILE_PATH, adata=None)
    """
    global GENE_COLOR_MAP, CELL_TYPE_COLORS

    if genes:
        make_color_map(genes)

    # get cell types from markers
    if os.path.exists(markers_file_path):
        cells_df = pd.read_csv(markers_file_path)
        cell_types = list(cells_df["source"].unique())
    else:
        cell_types = []

    # add from AnnData obs
    if adata is not None and "predicted_cell_type" in adata.obs:
        obs_types = list(adata.obs["predicted_cell_type"].cat.categories)
        cell_types = sorted(set(cell_types + obs_types))

    CELL_TYPE_COLORS = make_color_map(cell_types)

    # add defaults
    CELL_TYPE_COLORS["None"] = "#d3d3d3"
    CELL_TYPE_COLORS["BORIS marker"] = "#dd6c0f"

    for g in ["CTCF", "CTCFL", "Ctcfl", "Ctcflos"]:
        if g in GENE_COLOR_MAP:
            GENE_COLOR_MAP[f"{g}_pos"] = GENE_COLOR_MAP[g]   # positive = same as base
            GENE_COLOR_MAP[f"{g}_neg"] = "#d3d3d3"           # negative always grey

    return GENE_COLOR_MAP, CELL_TYPE_COLORS

def get_raw_gene_counts(adata, gene):
    # Priority is raw counts
    if adata.raw is not None:
        mat, genes = adata.raw.X, adata.var_names

    # Check for raw counts layer
    elif "counts" in adata.layers:
        mat, genes = adata.layers["counts"], adata.var_names

    # fallback to X if integer-like
    else:
        X = adata.X
        is_int_like = np.allclose(
            X[:100].toarray() if sp.issparse(X) else X[:100],
            (X[:100].toarray() if sp.issparse(X) else X[:100]).astype(int)
        )
        if not is_int_like:
            raise ValueError("No raw counts found (raw is None, no counts layer, X not integer).")
        mat, genes = X, adata.var_names

    if gene not in genes:
        raise KeyError(f"{gene} not found in gene list.")

    idx = list(genes).index(gene)
    col = mat[:, idx]
    return col.toarray().ravel() if sp.issparse(col) else np.ravel(col)

def get_cell_type(adata, markers_file_path=MARKERS_FILE_PATH, sample="sample"):
    """
    Compute cell-type enrichment (ULM) and annotate:
      - adata.obs['predicted_cell_type'] (per-cell max score)
      - adata.obs['cell_type'] (Leiden cluster → top cell type, fallback if needed)

      adata = get_cell_type(adata, markers_file_path=MARKERS_FILE_PATH, sample="sample")
    """

    
    if os.path.exists(markers_file_path):
        markers_df = pd.read_csv(markers_file_path, index_col=0)
        if not set(["source", "target"]).issubset(markers_df.columns):
            raise ValueError("markers_file must have columns: ['source','target']")
        markers_df = markers_df[["source", "target"]]
    else:
        
        try:
            markers = dc.get_resource(name="PanglaoDB", organism=organism, license="academic")
        except AttributeError:
            markers = dc.op.resource("PanglaoDB", organism=organism, license="academic")

        markers = markers[
            markers["human"].astype(bool)
            & markers["canonical_marker"].astype(bool)
            & (markers["human_sensitivity"].astype(float) > 0.5)
        ].copy()
        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
        markers_df = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})[["source", "target"]]
        markers_df.to_csv(markers_file_path)

    # normalize names
    adata.var_names = adata.var_names.astype(str).str.upper()
    net = markers_df.drop_duplicates(subset=['source', 'target']).copy()
    net['target'] = net['target'].astype(str).str.upper()

    # remove any duplicate targets after case-normalization
    net = net.drop_duplicates(subset=['target']).reset_index(drop=True)

    # overlap check
    overlap = set(adata.var_names) & set(net['target'])
    print(f"[{sample}] Overlapping genes: {len(overlap)}")

    if len(overlap) == 0:
        print(f"[{sample}] No overlap between adata.var_names and marker targets. Skipping ULM.")
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
        return adata

    # run ULM
    try:

        dc.mt.ulm(data=adata, net=net, tmin=0)
    except AssertionError as e:
        print(f"[{sample}] ULM failed: {e}. Using Leiden clusters as fallback.")
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
        return adata

    score = dc.pp.get_obsm(adata, key="score_ulm")
    score_df = score.to_df() if hasattr(score, "to_df") else score
    adata.obs["predicted_cell_type"] = score_df.idxmax(axis=1)

    # rank top cell types per cluster
    try:
        df = dc.tl.rankby_group(
            adata=score,
            groupby="leiden",
            reference="rest",
            method="t-test_overestim_var"
        )
    except Exception as e:
        print(f"[{sample}] rankby_group failed: {e}")
        df = None

    if df is None or df.empty:
        if "leiden" in adata.obs:
            adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
        else:
            print(f"[{sample}] No Leiden clusters found. Assigning all cells to 'unknown'.")
            adata.obs["cell_type"] = pd.Categorical(["unknown"] * adata.n_obs)

    else:
        df = df[df["stat"] > 0]
        dict_ann = (df.groupby("group").head(1)
                      .set_index("group")["name"].to_dict())
        s = adata.obs["leiden"].astype(str)
        adata.obs["cell_type"] = pd.Categorical(s.map(dict_ann).fillna(s))
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    return adata




