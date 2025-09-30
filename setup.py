# setup.py

#### Standard libraries ####
import os
import sys
import gzip
import shutil

#### Scientific stack ####
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.io import mmread
import random

#### Plotting ####
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorcet as cc

#### Domain-specific ####
import decoupler as dc

#### Custom ####
from config import *
#### Environmet Setup ####

# ---- globals ----
GENE_COLOR_MAP = {"None": "#d3d3d3"}
CELL_TYPE_COLORS = {}
GENE_COLOR_INDEX = 0


def make_color_map(items, fallback="#d3d3d3", palette="glasbey_light"):
    """
    Assign shuffled colors from a palette to items.
    Global state: GENE_COLOR_MAP, GENE_COLOR_INDEX
    """
    global GENE_COLOR_INDEX

    if palette in cc.palette:
        all_colors = list(cc.palette[palette])
    else:
        all_colors = list(getattr(cc, palette))

    random.seed(42)  # reproducible
    random.shuffle(all_colors)

    for g in items:
        if g not in GENE_COLOR_MAP:
            GENE_COLOR_MAP[g] = all_colors[GENE_COLOR_INDEX % len(all_colors)]
            GENE_COLOR_INDEX += 1

    # always enforce defaults
    GENE_COLOR_MAP["None"] = "#d3d3d3"
    GENE_COLOR_MAP["CTCFL"] = "#04531c"
    GENE_COLOR_MAP["CTCF"] = "#4b11d3"

    return {g: GENE_COLOR_MAP[g] for g in items}


def init_color_maps(genes=None, markers_file_path=MARKERS_FILE_PATH, adata=None):
    """
    Initialize gene and cell-type color maps.
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

    for g in ["CTCF", "CTCFL"]:
        if g in GENE_COLOR_MAP:
            GENE_COLOR_MAP[f"{g}_pos"] = GENE_COLOR_MAP[g]   # positive = same as base
            GENE_COLOR_MAP[f"{g}_neg"] = "#d3d3d3"           # negative always grey

    return GENE_COLOR_MAP, CELL_TYPE_COLORS

# Expand gene map to include _pos/_neg keys


#### Project Setup ####

#### Project Setup ####

def directory_setup(projects=None):
    """
    Create required directory structure for one or more projects.

    Usage:
        python main.py setup [project_name ...]
        If no project_name given, all projects from PROJECTS in config are used.
    """
    if projects is None:
        projects = list(PROJECTS.keys())

    # references directory
    if not os.path.exists(REFERENCE_DIR):
        os.makedirs(REFERENCE_DIR)
        print(f"[NEW] {REFERENCE_DIR}")
    else:
        print(f"[SKIP] {REFERENCE_DIR}")

    for project in projects:
        if project not in PROJECTS:
            print(f"[WARN] Unknown project: {project}, skipping")
            continue

        project_dir = os.path.join(BASE_DIR, project)

        new_dirs = [
            os.path.join(project_dir, RAW_DATA_DIR),
            os.path.join(project_dir, WORKING_DATA_DIR),
            os.path.join(project_dir, FIGURES_DIR),
            os.path.join(project_dir, TABLES_DIR),
        ]

        for d in new_dirs:
            if not os.path.exists(d):
                os.makedirs(d)
                print(f"[NEW] {d}")
            else:
                print(f"[SKIP] {d}")


def set_active_project(project_name):
    """
    Activate a project: return its directories, raw_file_type, and markers_file_path.
    Falls back to global MARKERS_FILE_PATH if not provided for the project.
    """
    if project_name not in PROJECTS:
        raise ValueError(f"Unknown project: {project_name}. Must be one of {list(PROJECTS.keys())}")

    raw_data_dir = os.path.join(BASE_DIR, project_name, RAW_DATA_DIR)
    working_adata_dir = os.path.join(BASE_DIR, project_name, WORKING_DATA_DIR)
    figures_dir = os.path.join(BASE_DIR, project_name, FIGURES_DIR)
    tables_dir = os.path.join(BASE_DIR, project_name, TABLES_DIR)

    raw_file_type = PROJECTS[project_name]["raw_file_type"]
    markers_file_path = PROJECTS[project_name].get("markers_file_path", MARKERS_FILE_PATH)

    return raw_data_dir, working_adata_dir, figures_dir, tables_dir, raw_file_type, markers_file_path



#### Raw Data Setup ####

def detect_compressed_files(raw_data_dir):
    """
    The function will detect and decompress .gz files in the given directory.

    Usage: detect_compressed_files(RAW_DATA_DIR)
    """
    for fname in os.listdir(raw_data_dir):
        full = os.path.join(raw_data_dir, fname)
        if fname.endswith(".gz"):
            with gzip.open(full, "rb") as f_in:
                with open(full[:-3], "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(full)


def add_list_to_config(list_name, list_values, config_path="config.py"):

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
    """
    Normalize gene annotations and set var_names to HGNC symbols.
    Preserves ensembl_id and hgnc_symbol, freezes counts into raw.
    """
    # keep hgnc_symbol + ensembl_id
    adata.var["hgnc_symbol"] = adata.var["gene_id"].astype(str)
    adata.var = adata.var[["ensembl_id", "hgnc_symbol"]]
    adata.var_names = adata.var["hgnc_symbol"]
    adata.var_names_make_unique()
    adata.var.index.name = None

    # freeze counts
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()

    print(adata)
    print("var head:\n", adata.var.head())

    # optional: adopt 'feature_name'
    if "feature_name" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names.astype(str)
        adata.var_names = adata.var["feature_name"].astype(str)
        adata.var.drop(columns=["feature_name"], inplace=True)

    # tidy indices
    adata.var.index.name = None
    if adata.raw is not None and adata.raw.var is not None:
        adata.raw.var.index.name = None

    return adata


def make_markers_df(markers_file_path=MARKERS_FILE_PATH):
    """
    Load or build a markers dataframe.
    Returns df with columns: source (cell type), target (gene).
    """
    if os.path.exists(markers_file_path):
        markers_df = pd.read_csv(markers_file_path, index_col=0)
    else:
        markers = dc.op.resource("PanglaoDB", organism="human")
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
def get_raw_gene_counts(adata, gene):
    # --- case 1: prefer raw ---
    if adata.raw is not None:
        mat, genes = adata.raw.X, adata.raw.var_names

    # --- case 2: counts layer ---
    elif "counts" in adata.layers:
        mat, genes = adata.layers["counts"], adata.var_names

    # --- case 3: fallback to X if integer-like ---
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


def get_raw_var_names(adata):
    """
    Return the gene names from the best available raw-count source.
    Priority:
      1. X if it looks integer-like
      2. counts layer
      3. adata.raw
    """
    # --- case 1: X looks like raw counts (integer-like) ---
    X = adata.X
    if sp.issparse(X):
        arr = X[:100].toarray()
    else:
        arr = np.asarray(X[:100])
    is_int_like = np.allclose(arr, arr.astype(int))

    if is_int_like:
        return adata.var_names

    # --- case 2: counts layer exists ---
    elif "counts" in adata.layers:
        return adata.var_names

    # --- case 3: fall back to raw ---
    elif adata.raw is not None:
        return adata.raw.var_names

    else:
        raise ValueError("No raw counts source found (X not integer, no counts layer, raw is None).")


def get_cell_type(adata, markers_file_path=MARKERS_FILE_PATH, sample="sample"):
    """
    Compute cell-type enrichment (ULM) and annotate:
      - adata.obs['predicted_cell_type'] (per-cell max score)
      - adata.obs['cell_type'] (Leiden cluster → top cell type, fallback if needed)
    """

    # ---- load or fetch marker network ----
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
        markers_df = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})[["source", "target"]]
        markers_df.to_csv(markers_file_path)

    # ---- normalize names ----
    adata.var_names = adata.var_names.astype(str).str.upper()
    net = markers_df.drop_duplicates(subset=['source', 'target']).copy()
    net['target'] = net['target'].astype(str).str.upper()

    # remove any duplicate targets after case-normalization
    net = net.drop_duplicates(subset=['target']).reset_index(drop=True)

    # ---- overlap check ----
    overlap = set(adata.var_names) & set(net['target'])
    print(f"[{sample}] Overlapping genes: {len(overlap)}")

    if len(overlap) == 0:
        print(f"[{sample}] No overlap between adata.var_names and marker targets. Skipping ULM.")
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
        return adata

    # ---- run ULM ----
    try:
        dc.mt.ulm(data=adata, net=net, tmin=2)
    except AssertionError as e:
        print(f"[{sample}] ULM failed: {e}. Using Leiden clusters as fallback.")
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str)
        return adata

    score = dc.pp.get_obsm(adata, key="score_ulm")
    score_df = score.to_df() if hasattr(score, "to_df") else score
    adata.obs["predicted_cell_type"] = score_df.idxmax(axis=1)

    # ---- rank top cell types per cluster ----
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



