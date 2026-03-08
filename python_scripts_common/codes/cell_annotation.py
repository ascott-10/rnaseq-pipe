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

def get_cell_type(sample, adata, markers_file_path, cell_type_label_col="predicted_cell_type"):
    

    if os.path.exists(markers_file_path):
        markers_df = pd.read_csv(markers_file_path, index_col=0)
        if not set(["source", "target"]).issubset(markers_df.columns):
            raise ValueError("markers_file must have columns: ['source','target']")
        markers_df = markers_df[["source", "target"]]
    else:
        try:
            markers = dc.get_resource(name="PanglaoDB", organism="human", license="academic")
        except AttributeError:
            markers = dc.op.resource("PanglaoDB", organism="human", license="academic")

        markers = markers[
            markers["human"].astype(bool)
            & markers["canonical_marker"].astype(bool)
            & (markers["human_sensitivity"].astype(float) > 0.5)
        ].copy()

        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

        markers_df = markers.rename(
            columns={"cell_type": "source", "genesymbol": "target"}
        )[["source", "target"]]

        markers_df.to_csv(markers_file_path)

    # normalize gene names
    adata.var_names = adata.var_names.astype(str).str.upper()

    net = markers_df.drop_duplicates(subset=["source", "target"]).copy()
    net["target"] = net["target"].astype(str).str.upper()
    net = net.reset_index(drop=True)

    # check overlap
    overlap = set(adata.var_names) & set(net["target"])
    print(f"[{sample}] Overlapping genes: {len(overlap)}")

    if len(overlap) == 0:
        print(f"[{sample}] No overlap between adata.var_names and marker targets. Skipping ULM.")
        adata.obs[cell_type_label_col] = adata.obs["leiden"].astype(str)
        return adata

    # run ULM
    try:
        dc.mt.ulm(data=adata, net=net, tmin=0)
    except AssertionError as e:
        print(f"[{sample}] ULM failed: {e}. Using Leiden clusters as fallback.")
        adata.obs[cell_type_label_col] = adata.obs["leiden"].astype(str)
        return adata

    # get ULM scores
    score = dc.pp.get_obsm(adata, key="score_ulm")
    score_df = score.to_df() if hasattr(score, "to_df") else score

    # per-cell predicted cell type
    adata.obs["ulm_score_cell_type"] = score_df.idxmax(axis=1)

    # cluster-level ranking
    try:
        df = dc.tl.rankby_group( adata=score, groupby="leiden", reference="rest",method="t-test_overestim_var")
    except Exception as e:
        print(f"[{sample}] rankby_group failed: {e}")
        df = None

    if df is None or df.empty:

        if "leiden" in adata.obs:
            adata.obs[cell_type_label_col] = adata.obs["leiden"].astype(str)
        else:
            print(f"[{sample}] No Leiden clusters found. Assigning all cells to 'unknown'.")
            adata.obs[cell_type_label_col] = pd.Categorical(["unknown"] * adata.n_obs)

    else:
        df = df[df["stat"] > 0]
        dict_ann = (df.groupby("group").head(1).set_index("group")["name"].to_dict())
        s = adata.obs["leiden"].astype(str)
        adata.obs[cell_type_label_col] = pd.Categorical(s.map(dict_ann).fillna(s))

    adata.obs[cell_type_label_col] = adata.obs[cell_type_label_col].astype("category")

    return adata


def annotate_cells(adata, cell_annotation_file_path, cell_type_label_col="predicted_cell_type"):

    cellmarkers_annotation_df = pd.read_csv(cell_annotation_file_path)

    hierarchy_cols = ["cell_name","cell_type","tissue_class","tissue_type", "cancer_type"]

    def most_common(x):
        return x.mode().iloc[0] if not x.mode().empty else None

    annotation_df = (
        cellmarkers_annotation_df[hierarchy_cols]
        .groupby("cell_name", as_index=False)
        .agg(most_common)
    )



    annotation_df[cell_type_label_col] = annotation_df["cell_name"]

    # set lookup index
    annotation_df = annotation_df.set_index(cell_type_label_col)

    for col in annotation_df.columns:

        if col == cell_type_label_col:
            continue

        adata.obs[col] = adata.obs[cell_type_label_col].map(annotation_df[col])

    return adata

def add_gene_binary_columns(adata, genes, threshold=0.1):

    if adata.raw is None:
        raise ValueError("adata.raw is None")

    for g in genes:

        if g not in adata.raw.var_names:
            continue

        expr = adata.raw[:, g].X

        if sp.issparse(expr):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()

        pos_mask = expr > threshold

        adata.obs[f"{g}_binary"] = pd.Categorical(
            np.where(pos_mask, f"{g}_pos", f"{g}_neg"),
            categories=[f"{g}_neg", f"{g}_pos"],
        )

    return adata

def set_obs_colors(adata, col, palette, cell_type_colors, gene_color_map):

    import colorcet as cc
    

    # convert palette name
    if isinstance(palette, str):
        palette = list(cc.palette[palette])

    adata.obs[col] = adata.obs[col].astype("category")

    # stable category order
    categories = sorted(adata.obs[col].cat.categories)

    assigned = []

    for i, cat in enumerate(categories):

        # predicted cell types
        if col == "predicted_cell_type":

            if cat not in cell_type_colors:
                cell_type_colors[cat] = palette[i % len(palette)]

            assigned.append(cell_type_colors[cat])

        # gene binary
        elif col.endswith("_binary"):

            gene = col.replace("_binary", "")

            if cat == f"{gene}_pos":

                if gene not in gene_color_map:
                    gene_color_map[gene] = palette[i % len(palette)]

                assigned.append(gene_color_map[gene])

            else:
                assigned.append("#d3d3d3")

        else:
            assigned.append(palette[i % len(palette)])

    adata.uns[f"{col}_colors"] = assigned