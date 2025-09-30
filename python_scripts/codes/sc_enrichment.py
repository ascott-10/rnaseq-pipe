# python_scripts/codes/sc_enrichment.py

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

def tf_scoring(adata, sample, figures_dir, tf_list = ["SMAD1", "SMAD2", "SMAD3", "SMAD4"], plot_yes = False):
    os.makedirs(os.path.join(figures_dir, "collectri"), exist_ok=True)
    
    
    """CollecTRI network
CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources [MullerDTV+23]. This collection provides an increased coverage of transcription factors and a superior performance in identifying perturbed TFs compared to other literature based GRNs such as DoRothEA [GAHI+19]. Similar to DoRothEA, interactions are weighted by their mode of regulation (activation or inhibition)."""

    # CollecTRI network
    collectri = dc.op.collectri(organism="human")

    # TF scores (persist under unique key)
    if "score_ulm_tf" not in adata.obsm:
        dc.mt.ulm(data=adata, net=collectri)  # writes to 'score_ulm'
        s = dc.pp.get_obsm(adata, key="score_ulm")
        tf_df = s.to_df() if hasattr(s, "to_df") else s
        adata.obsm["score_ulm_tf"] = tf_df
        if "score_ulm" in adata.obsm:
            del adata.obsm["score_ulm"]
    else:
        tf_scores = adata.obsm["score_ulm_tf"]
        tf_df = tf_scores.to_df() if hasattr(tf_scores, "to_df") else tf_scores

    # rank markers per cell type (work on a lightweight AnnData built from tf_df)
    sadata = anndata.AnnData(
        X=tf_df.values,
        obs=adata.obs[["cell_type"]].copy(),
        var=pd.DataFrame(index=tf_df.columns),
    )
    sadata.obs["cell_type"] = sadata.obs["cell_type"].astype("category")

    df = dc.tl.rankby_group(
        adata=sadata, groupby="cell_type", reference="rest", method="t-test_overestim_var"
    )
    df = df[df["stat"] > 0]

    n_markers = 3
    source_markers = (
        df.groupby("group")
          .head(n_markers)
          .drop_duplicates("name")
          .groupby("group")["name"]
          .apply(lambda x: list(x))
          .to_dict()
    )

    # per-cell top TF label + its score
    adata.obs["top_tf"] = tf_df.idxmax(axis=1).astype("category")
    adata.obs["top_tf_score"] = tf_df.max(axis=1).astype(float)

    if plot_yes:
        # UMAP of TF gene expression + cell_type
        colors = list(tf_list) + ["cell_type"]
        sc.pl.umap(adata, color=colors, frameon=False, show=False)
        plt.savefig(os.path.join(figures_dir, "collectri", f"{sample}_tfs_expression_umap.png"),
                    bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.close()

        # Matrixplot of top TF markers per cell type (on score space)
        sc.pl.matrixplot(
            adata=sadata,
            var_names=source_markers,
            groupby="cell_type",
            dendrogram=True,
            standard_scale="var",
            colorbar_title="Z-scaled scores",
            cmap="RdBu_r",
            show=False,
        )
        plt.savefig(os.path.join(figures_dir, "collectri", f"{sample}_tfs_celltype_matrix.png"),
                    bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.close()

        # Network plot (mean expression + mean TF scores)
        # pick the strongest TFs to show (by mean TF score)
        # pick top 2 TFs by mean activity
        top_tfs = adata.obsm["score_ulm_tf"].mean().nlargest(2).index.tolist()
        gex   = adata.to_df().mean(0).to_frame().T
        score = adata.obsm["score_ulm_tf"].mean(0).to_frame().T

        dc.pl.network(
            data=gex, score=score, net=dc.op.collectri(organism="human"),
            sources=top_tfs, targets=10, size_node=7, figsize=(8,5), s_cmap="Reds"
        )
        plt.setp(plt.gca().texts, fontsize=6)  # shrink labels (one-liner)
        plt.savefig(os.path.join(figures_dir, "collectri", f"{sample}_cells_tfs_networks.png"),
                    bbox_inches="tight", dpi=300)
        plt.close()

        # compress top TFs to 12; rest -> "Other"
        # show only the K most represented TFs in the legend
        K = 12
        s = adata.obs["top_tf"].astype(str)
        keep = s.value_counts().index[:K].tolist()
        s = s.where(s.isin(keep), "Other")
        cats = keep + (["Other"] if "Other" in s.unique() else [])
        adata.obs["top_tf_top"] = pd.Categorical(s, categories=cats)

        # plot
        sc.pl.umap(
            adata,
            color=["top_tf_top", "cell_type"],
            legend_loc="right margin",   # legends on the right for both panels
            frameon=False,
            show=False,
        )
        plt.gcf().set_size_inches(10, 4)
        plt.savefig(os.path.join(figures_dir, "collectri", f"{sample}_umap_top_tf_with_legends.png"),
                    bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.close()

    return adata

def pathway_scoring(adata, sample, figures_dir, path_list = ["NFkB", "TGFb", "Hypoxia"], plot_yes = False):
    os.makedirs(os.path.join(figures_dir, "pathways"), exist_ok=True)
    os.makedirs(os.path.join(figures_dir, "CTCFL"), exist_ok=True)
    
    """PROGENy is a comprehensive resource that provides a curated collection of pathways and their target genes, along with weights for each interaction [SKK+18]. This collection is designed to infer pathway activities from gene expression data, allowing researchers to gain insights into the underlying signaling pathways that drive cellular behavior. PROGENy covers a wide range of pathways, including those involved in cancer, immune response, and development."""


    progeny = dc.op.progeny(organism="human")  # ["source","target","weight","padj"]

    # compute once and store under a unique key
    if "score_ulm_progeny" not in adata.obsm:
        dc.mt.ulm(data=adata, net=progeny)             # writes to 'score_ulm'
        s = dc.pp.get_obsm(adata, key="score_ulm")
        pw_df = s.to_df() if hasattr(s, "to_df") else s
        adata.obsm["score_ulm_progeny"] = pw_df
        if "score_ulm" in adata.obsm:   # avoid future collisions
            del adata.obsm["score_ulm"]
    else:
        tmp = adata.obsm["score_ulm_progeny"]
        pw_df = tmp.to_df() if hasattr(tmp, "to_df") else tmp

    # copy selected pathway scores into obs for UMAP coloring
    for p in path_list:
        if p in pw_df.columns:
            adata.obs[f"{p}_score"] = pw_df[p]

    if plot_yes:
        # UMAP (plot on adata, which has X_umap)
        colors = list(path_list) + ["cell_type"]
        sc.pl.umap(
            adata,
            color=[f"{c}_score" if c in pw_df.columns else c for c in colors],
            cmap="RdBu_r",
            vcenter=0,
            frameon=False,
            show=False,
        )
        plt.savefig(os.path.join(figures_dir, "pathways", f"{sample}_pathways_umap.png"),
                    bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.close()

        # matrixplot on a lightweight AnnData built from the score matrix
        sadata = anndata.AnnData(
            X=pw_df.values,
            obs=adata.obs[["cell_type"]].copy(),
            var=pd.DataFrame(index=pw_df.columns),
        )
        sadata.obs["cell_type"] = sadata.obs["cell_type"].astype("category")

        # show ALL pathways
        sc.pl.matrixplot(
            adata=sadata,
            var_names=sadata.var_names,   # <-- all pathways
            groupby="cell_type",
            dendrogram=True,
            standard_scale="var",
            colorbar_title="Z-scaled scores",
            cmap="RdBu_r",
            show=False,
        )
        plt.savefig(os.path.join(figures_dir, "pathways", f"{sample}_pathways_matrix.png"),
                    bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.close()



    return adata

def hallmark_genesets(adata, sample, figures_dir, plot_yes = False):
    os.makedirs(os.path.join(figures_dir, "genesets"), exist_ok=True)
    """Hallmark gene sets are curated collections of genes that represent specific, well-defined biological states or processes. They are part of MSigDB and were developed to reduce redundancy and improve interpretability compared to older, more overlapping gene set collections [LSP+11].

    A total of 50 gene sets are provided, designed to be non-redundant, concise, and biologically coherent."""


    hallmark = dc.op.hallmark(organism="human")  # ["source","target"]

    # compute once and store under a unique key
    if "score_ulm_hallmark" not in adata.obsm:
        dc.mt.ulm(data=adata, net=hallmark)               # writes to 'score_ulm'
        s = dc.pp.get_obsm(adata, key="score_ulm")
        hk_df = s.to_df() if hasattr(s, "to_df") else s
        adata.obsm["score_ulm_hallmark"] = hk_df
        if "score_ulm" in adata.obsm:                     # avoid collisions later
            del adata.obsm["score_ulm"]
    else:
        tmp = adata.obsm["score_ulm_hallmark"]
        hk_df = tmp.to_df() if hasattr(tmp, "to_df") else tmp

    # rank markers per cell type on score space
    sadata = anndata.AnnData(
        X=hk_df.values,
        obs=adata.obs[["cell_type"]].copy(),
        var=pd.DataFrame(index=hk_df.columns),
    )
    sadata.obs["cell_type"] = sadata.obs["cell_type"].astype("category")

    df = dc.tl.rankby_group(adata=sadata, groupby="cell_type",
                            reference="rest", method="t-test_overestim_var")
    df = df[df["stat"] > 0]

    n_markers = 3
    source_markers = (
        df.groupby("group")
          .head(n_markers)
          .drop_duplicates("name")
          .groupby("group")["name"]
          .apply(list)
          .to_dict()
    )

    if plot_yes:
        sc.pl.matrixplot(
            adata=sadata,
            var_names=source_markers,
            groupby="cell_type",
            dendrogram=True,
            standard_scale="var",
            colorbar_title="Z-scaled scores",
            cmap="RdBu_r",
            show=False,
        )
        plt.savefig(os.path.join(figures_dir, "genesets", f"{sample}_hallmark_genesets_matrix.png"),
                    bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.close()

    return adata

def ctcfl_pathway_matrices(adata, sample, figures_dir):
    os.makedirs(os.path.join(figures_dir, "CTCFL"), exist_ok=True)

    # PROGENy scores (cells x pathways)
    if "score_ulm_progeny" in adata.obsm:
        pw_df = adata.obsm["score_ulm_progeny"]
        pw_df = pw_df.to_df() if hasattr(pw_df, "to_df") else pw_df
    else:
        progeny = dc.op.progeny(organism="human")
        dc.mt.ulm(data=adata, net=progeny)
        s = dc.pp.get_obsm(adata, key="score_ulm")
        pw_df = s.to_df() if hasattr(s, "to_df") else s
        adata.obsm["score_ulm_progeny"] = pw_df
        if "score_ulm" in adata.obsm: del adata.obsm["score_ulm"]

    # split by labels in obs['CTCFL_expr']
    lab = adata.obs["CTCFL_expr"].astype(str).str.lower()
    pos = lab.isin(["CTCFL_pos", "pos", "positive", "CTCFL+"])
    neg = lab.isin(["CTCFL_neg", "neg", "negative", "CTCFL-"])

    # score-space AnnData
    s_pos = anndata.AnnData(X=pw_df.loc[pos].values,
                            obs=adata.obs.loc[pos, ["cell_type"]].copy(),
                            var=pd.DataFrame(index=pw_df.columns))
    s_neg = anndata.AnnData(X=pw_df.loc[neg].values,
                            obs=adata.obs.loc[neg, ["cell_type"]].copy(),
                            var=pd.DataFrame(index=pw_df.columns))
    for s in (s_pos, s_neg):
        s.obs["cell_type"] = s.obs["cell_type"].astype("category")
        if hasattr(s.obs["cell_type"], "cat"):
            s.obs["cell_type"] = s.obs["cell_type"].cat.remove_unused_categories()

    if s_pos.n_obs == 0 or s_neg.n_obs == 0:
        print("[warn] one CTCFL group has 0 cells; skipping pathways matrixplot")
        return adata

    # save matrices (ALL pathways)
    pdir = os.path.join(figures_dir, "CTCFL")
    p_pos = os.path.join(pdir, f"{sample}_pathways_matrix_CTCFLpos.png")
    p_neg = os.path.join(pdir, f"{sample}_pathways_matrix_CTCFLneg.png")

    sc.tl.dendrogram(s_pos, groupby="cell_type")
    sc.tl.dendrogram(s_neg, groupby="cell_type")

    sc.pl.matrixplot(s_pos, var_names=s_pos.var_names, groupby="cell_type",
                     dendrogram=True, standard_scale="var", colorbar_title="Z-scaled scores",
                     cmap="RdBu_r", show=False)
    plt.savefig(p_pos, bbox_inches="tight", pad_inches=0.1, dpi=300); plt.close()

    sc.pl.matrixplot(s_neg, var_names=s_neg.var_names, groupby="cell_type",
                     dendrogram=True, standard_scale="var", colorbar_title="Z-scaled scores",
                     cmap="RdBu_r", show=False)
    plt.savefig(p_neg, bbox_inches="tight", pad_inches=0.1, dpi=300); plt.close()

    

    # side-by-side composite (PROGENy)
    import matplotlib.image as mpimg
    img1, img2 = mpimg.imread(p_pos), mpimg.imread(p_neg)
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    ax[0].imshow(img1); ax[0].axis("off"); ax[0].set_title("CTCFL_pos (PROGENy)")
    ax[1].imshow(img2); ax[1].axis("off"); ax[1].set_title("CTCFL_neg (PROGENy)")
    plt.tight_layout()
    plt.savefig(os.path.join(pdir, f"{sample}_CTCFLpos_vs_neg_pathways_matrix.png"),
                bbox_inches="tight", dpi=300)
    plt.close()

    return adata


def ctcfl_geneset_matrices(adata, sample, figures_dir):
    os.makedirs(os.path.join(figures_dir, "CTCFL"), exist_ok=True)

    # Hallmark scores (cells x genesets)
    if "score_ulm_hallmark" in adata.obsm:
        hk_df = adata.obsm["score_ulm_hallmark"]
        hk_df = hk_df.to_df() if hasattr(hk_df, "to_df") else hk_df
    else:
        hallmark = dc.op.hallmark(organism="human")
        dc.mt.ulm(data=adata, net=hallmark)
        s = dc.pp.get_obsm(adata, key="score_ulm")
        hk_df = s.to_df() if hasattr(s, "to_df") else s
        adata.obsm["score_ulm_hallmark"] = hk_df
        if "score_ulm" in adata.obsm: del adata.obsm["score_ulm"]

    # split by labels in obs['CTCFL_expr']
    lab = adata.obs["CTCFL_expr"].astype(str).str.lower()
    pos = lab.isin(["CTCFL_pos", "pos", "positive", "CTCFL+"])
    neg = lab.isin(["CTCFL_neg", "neg", "negative", "CTCFL-"])

    # score-space AnnData
    s_pos = anndata.AnnData(X=hk_df.loc[pos].values,
                            obs=adata.obs.loc[pos, ["cell_type"]].copy(),
                            var=pd.DataFrame(index=hk_df.columns))
    s_neg = anndata.AnnData(X=hk_df.loc[neg].values,
                            obs=adata.obs.loc[neg, ["cell_type"]].copy(),
                            var=pd.DataFrame(index=hk_df.columns))
    for s in (s_pos, s_neg):
        s.obs["cell_type"] = s.obs["cell_type"].astype("category")
        if hasattr(s.obs["cell_type"], "cat"):
            s.obs["cell_type"] = s.obs["cell_type"].cat.remove_unused_categories()

    if s_pos.n_obs == 0 or s_neg.n_obs == 0:
        print("[warn] one CTCFL group has 0 cells; skipping genesets matrixplot")
        return adata

    # save matrices (ALL genesets)
    pdir = os.path.join(figures_dir, "CTCFL")
    p_pos = os.path.join(pdir, f"{sample}_genesets_matrix_CTCFLpos.png")
    p_neg = os.path.join(pdir, f"{sample}_genesets_matrix_CTCFLneg.png")

    sc.tl.dendrogram(s_pos, groupby="cell_type")
    sc.tl.dendrogram(s_neg, groupby="cell_type")

    sc.pl.matrixplot(s_pos, var_names=s_pos.var_names, groupby="cell_type",
                     dendrogram=True, standard_scale="var", colorbar_title="Z-scaled scores",
                     cmap="RdBu_r", show=False)
    plt.savefig(p_pos, bbox_inches="tight", pad_inches=0.1, dpi=300); plt.close()

    sc.pl.matrixplot(s_neg, var_names=s_neg.var_names, groupby="cell_type",
                     dendrogram=True, standard_scale="var", colorbar_title="Z-scaled scores",
                     cmap="RdBu_r", show=False)
    plt.savefig(p_neg, bbox_inches="tight", pad_inches=0.1, dpi=300); plt.close()

    # side-by-side composite (Hallmark)
    import matplotlib.image as mpimg
    img1, img2 = mpimg.imread(p_pos), mpimg.imread(p_neg)
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    ax[0].imshow(img1); ax[0].axis("off"); ax[0].set_title("CTCFL_pos (Hallmark)")
    ax[1].imshow(img2); ax[1].axis("off"); ax[1].set_title("CTCFL_neg (Hallmark)")
    plt.tight_layout()
    plt.savefig(os.path.join(pdir, f"{sample}_CTCFLpos_vs_neg_genesets_matrix.png"),
                bbox_inches="tight", dpi=300)
    plt.close()

    return adata
#End of sc_enrichment.py