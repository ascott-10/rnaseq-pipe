# projects/E-MTAB-8559/python_scripts/codes/de_expression_py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc
import gseapy as gp

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread
from scipy.stats import spearmanr

import matplotlib.pyplot as plt

#### Custom ####

from config import *
from codes.cell_annotation import make_universal_gene_list, init_color_maps

### Make palettes ### 
CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)
GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CLEANED_GENES,markers_file_path=MARKERS_FILE_PATH, adata=None)

### Functions ###

def rank_genes(sample, adata_init_path, adata_ranked_path, gene="CTCFL"):
    """adata_ranked = rank_genes(sample, adata_init_path, adata_ranked_path, gene="CTCFL")"""
    adata = sc.read_h5ad(adata_init_path)

    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None:
        adata.X = adata.raw.X.copy()
        adata.var = adata.raw.var.copy()

    if f"{gene}_binary" not in adata.obs.columns:
        gene_idx = np.where(adata.var_names == gene)[0][0]
        expr = adata.raw.X[:, gene_idx].toarray().ravel()
        adata.obs[f"{gene}_binary"] = np.where(expr > 0, "pos", "neg")
        adata.obs[f"{gene}_binary"] = adata.obs[f"{gene}_binary"].astype("category")

    adata.obs[f"{gene}_binary"] = (
        adata.obs[f"{gene}_binary"]
        .cat.remove_unused_categories()
    )

    counts = adata.obs[f"{gene}_binary"].value_counts()

    if "pos" not in counts or "neg" not in counts:
        if "rank_genes_groups" in adata.uns:
            del adata.uns["rank_genes_groups"]
        return adata

    if counts["pos"] < 2 or counts["neg"] < 2:
        if "rank_genes_groups" in adata.uns:
            del adata.uns["rank_genes_groups"]
        return adata

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    sc.tl.rank_genes_groups(
    adata,
    groupby="CTCFL_binary",
    method="wilcoxon"
)
    adata.write(adata_ranked_path)
      
    return adata


def plot_ranked_genes(adata_ranked_path, sample, fig_dir = FIGURES_DIR, tables_dir = TABLES_DIR, group_key="CTCFL_binary", group_name="pos", groupby="label", fdr=0.01, log_fold_change=1):
    
    table_save_path = os.path.join(tables_dir, "DE_analysis")
    fig_save_path = os.path.join(fig_dir, "DE_analysis")
    os.makedirs(fig_save_path, exist_ok=True)
    os.makedirs(table_save_path, exist_ok=True)
    sc.settings.figdir = fig_save_path

    adata = sc.read_h5ad(adata_ranked_path)

    

    if group_key not in adata.obs:
        print(f"{group_key} not in adata.obs")
        return

    cats = adata.obs[group_key].astype("category").cat.categories

    if group_name not in cats:
        print(f"{group_name} not in {group_key}")
        return

    if len(cats) < 2:
        print(f"len({group_key}) < 2")
        return

    if "rank_genes_groups" not in adata.uns:
        print("No rank_genes_groups found, skipping")
        return

    

    sc.pl.rank_genes_groups_dotplot(adata, 
                                groupby=group_key, 
                               
                                n_genes=20,
                                values_to_plot="logfoldchanges", 
                                show=True, 
                                save=f"{sample}_{group_name}_dotplot.png")

    res = sc.get.rank_genes_groups_df( adata, group=group_name)
    res.index = res["names"].values
    res = res[ (res["pvals_adj"] < fdr) & (res["logfoldchanges"] > log_fold_change)].sort_values(by="logfoldchanges")

    
    max_genes = 30
    markers = res.nlargest(max_genes, "logfoldchanges").index.tolist()
    print(f"Plotting {len(res)} genes...")

    sc.pl.heatmap(
    adata,
    markers,
    groupby="CTCFL_binary",
    swap_axes=True,
    show=True,
    save=f"{sample}_{group_name}_heatmap.png"
)

    res.to_csv(os.path.join(table_save_path, f"{sample}_rank_genes_{group_name}.csv"),index=False)

def correlate_with_gene(adata_init_path, sample, fig_dir=FIGURES_DIR, tables_dir=TABLES_DIR, n_top=30, gene="CTCFL", method="spearman"):

    table_save_path = os.path.join(tables_dir, "DE_analysis")
    fig_save_path = os.path.join(fig_dir, "DE_analysis")
    os.makedirs(fig_save_path, exist_ok=True)
    os.makedirs(table_save_path, exist_ok=True)
    sc.settings.figdir = fig_save_path

    adata = sc.read_h5ad(adata_init_path)

    if "counts" in adata.layers:
        X = adata.layers["counts"]
    elif adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    if gene not in adata.var_names:
        raise ValueError(f"{gene} not found in var_names")

    target_idx = np.where(adata.var_names == gene)[0][0]
    target_expr = X[:, target_idx].toarray().ravel()

    if np.all(target_expr == 0):
        print(f"{sample}: no expression of {gene}, skipping correlation")
        return


    corrs = []
    pvals = []

    for i, g in enumerate(adata.var_names):
        gene_expr = X[:, i].toarray().ravel()

        if method == "spearman":
            if np.std(gene_expr) == 0 or np.std(target_expr) == 0:
                r, p = np.nan, np.nan
            else:
                r, p = spearmanr(target_expr, gene_expr)
        else:
            r = np.corrcoef(target_expr, gene_expr)[0, 1]
            p = np.nan

        corrs.append(r)
        pvals.append(p)

    res = pd.DataFrame(
        {
            "gene": adata.var_names,
            "correlation": corrs,
            "pval": pvals,
        }
    ).sort_values("correlation", ascending=False)

    res.to_csv(
        os.path.join(table_save_path, f"{sample}_{gene}_correlations.csv"),
        index=False
    )

    top_genes = res["gene"].head(n_top).tolist()

    adata.obs["_all"] = "all"

    sc.pl.heatmap(
        adata,
        top_genes,
        groupby="_all",
        swap_axes=True,
        show=True,
        save=f"{sample}_{gene}_correlation_heatmap.png"
    )

    adata.obs.drop(columns="_all", inplace=True)


def run_cluster_signature(
    bdata,
    cluster_id,
    sample,
    outdir,
    groupby="CTCFL_binary",
    method="wilcoxon",
    gene_sets_gsea="KEGG_2016",
    gene_sets_go="GO_Biological_Process_2021"
):
    sample_dir = os.path.join(outdir, sample)
    os.makedirs(sample_dir, exist_ok=True)

    cluster_label = str(cluster_id)
    groups = bdata.obs[groupby].astype(str)

    if cluster_label not in groups.unique():
        print(f"[{sample}] Cluster {cluster_label} not found, skipping")
        return None

    n_cells = (groups == cluster_label).sum()
    if n_cells < 5:
        print(f"[{sample}] Cluster {cluster_label} has {n_cells} cells, skipping")
        return None

    sc.tl.rank_genes_groups(
        bdata,
        groupby=groupby,
        groups=[cluster_label],
        reference="rest",
        method=method
    )

    result = bdata.uns["rank_genes_groups"]

    degs = pd.DataFrame({
        "gene": result["names"][cluster_label],
        "score": result["scores"][cluster_label],
        "pval": result["pvals"][cluster_label],
        "pval_adj": result["pvals_adj"][cluster_label],
        "logfoldchange": result["logfoldchanges"][cluster_label],
    }).dropna()

    degs.to_csv(
    f"{sample_dir}/cluster{cluster_id}_DEGs.csv",
    index=False
    )

    pre_res = gp.prerank(
        rnk=degs[["gene", "logfoldchange"]]
            .set_index("gene")
            .squeeze(),
        gene_sets=gene_sets_gsea
    )

    degs_sig = degs[degs["pval_adj"] < 0.05]
    degs_up = degs_sig[degs_sig["logfoldchange"] > 0]
    degs_dw = degs_sig[degs_sig["logfoldchange"] < 0]

    enr_up = None
    enr_dw = None

    if len(degs_up) > 0:
        enr_up = gp.enrichr(
            gene_list=degs_up["gene"].tolist(),
            gene_sets=gene_sets_go,
            outdir=None
        )

    if len(degs_dw) > 0:
        enr_dw = gp.enrichr(
            gene_list=degs_dw["gene"].tolist(),
            gene_sets=gene_sets_go,
            outdir=None
        )

    go_tables = []

    if enr_up is not None:
        enr_up.res2d["Term"] = enr_up.res2d["Term"].str.split(" \\(GO").str[0]
        enr_up.res2d["UP_DW"] = "UP"
        go_tables.append(enr_up.res2d)

    if enr_dw is not None:
        enr_dw.res2d["Term"] = enr_dw.res2d["Term"].str.split(" \\(GO").str[0]
        enr_dw.res2d["UP_DW"] = "DOWN"
        go_tables.append(enr_dw.res2d)

    if len(go_tables) == 0:
        print(f"[{sample}][cluster {cluster_id}] No GO enrichment")
        return {
            "degs": degs,
            "gsea": pre_res.res2d,
            "go": None
        }

    enr_res = pd.concat(go_tables)
    
    enr_res.to_csv(
    f"{sample}_{cluster_id}_GO_UP_DOWN.csv",
        index=False
    )

    ax = gp.dotplot(
        enr_res,
        figsize=(5, 5),
        x="UP_DW",
        x_order=["UP", "DOWN"],
        title=f"Cluster {cluster_id} GO_BP",
        size=3,
        show_ring=False,
        show=False
    )

    ax.set_xlabel("Regulation direction")
    ax.set_ylabel("GO Biological Process")

    fig = ax.figure
    fig.subplots_adjust(right=0.75)

    fig.savefig(
    f"{sample}_{cluster_id}_GO_BP_dotplot.png",
    bbox_inches="tight",
    dpi=300)
    plt.close(fig)

    return {
        "degs": degs,
        "gsea": pre_res.res2d,
        "go": enr_res
    }
