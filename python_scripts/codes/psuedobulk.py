# python_scripts/codes/pseudobulk.py

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



# Project configuration (paths, constants, etc.)
from config import *
import gzip
import shutil
from codes.utils import make_markers_df, get_raw_gene_counts, get_cell_type, add_list_to_config

def map_metadata(adata,sample,adata_working_path,meta_columns=("tumor_stage", "chemotherapy_response"), reference_dir=REFERENCE_DIR):
    meta_path = os.path.join(reference_dir, "sample_metadata.csv")
    metadata_df = pd.read_csv(meta_path)
    metadata_df = metadata_df.set_index("gse_id")

    # handle both "T77" and "Tumor_T77"
    idx = metadata_df.index.astype(str)
    m = idx.str.contains(sample, regex=False); print(idx[m].tolist())
    key = idx[m][0]


    
    # ensure per-cell sample column exists (good practice)
    if "sample" not in adata.obs:
        adata.obs["sample"] = key

    # map each requested column (broadcast scalar to all cells)
    for col in meta_columns:
        if col in metadata_df.columns:
            adata.obs[col] = metadata_df.loc[key, col]
        else:
            print(f"[map_metadata] column '{col}' not in metadata CSV")

    # tidy types for plotting
    for col in ("tumor_stage", "chemotherapy_response"):
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")

    

    # tidy types
    for col in ("tumor_stage", "chemotherapy_response"):
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype("category")

   
        

    print(adata.obs.columns)
    print(f"Saving processed data to {adata_working_path}")
    adata.write(adata_working_path)
    return adata


import anndata
import os


def setup_pseudobulk(adata, sample, figures_dir, tables_dir):

    """The pseudo-bulk approach involves the following steps:

    1. Subsetting the donor groups of interest
    2. Extracting their raw integer counts
    3. Summing their counts per gene into a single profile if they pass quality control

    Then, DEA can be performed if there are at least two biological replicates per condition.

    Pseudobulking can easily be performed using the function decoupler.pp.pseudobulk().
    In this example, the counts are just summed, though other modes such as the mean or any custom aggregation
    function are available. For more information, refer to the mode argument.
    """

    os.makedirs(os.path.join(figures_dir, "pseudobulk"), exist_ok=True)
    os.makedirs(os.path.join(tables_dir, "pseudobulk"), exist_ok=True)
    print("Setting up pseudobulk analysis...")

    # --- run pseudobulk
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col="donor_id",                # 👈 changed from "sample" to "donor_id"
        groups_col=["CTCFL_expr", "cell_type"],  # still show cell type info for context
        layer="counts",
        mode="sum"
    )
    print("Initial pseudobulk profiles:", pdata.shape)

    # --- quick QC stats
    summary = (
        pdata.obs.groupby(["CTCFL_expr", "cell_type"])[["psbulk_counts", "psbulk_cells"]]
        .sum()
        .reset_index()
    )

    print(summary.head())   # quick check
    
    summary_save_path = os.path.join(tables_dir, "pseudobulk", f"{sample}_summary_pseudobulk.csv")
    summary_barplot_save_path = os.path.join(figures_dir, "pseudobulk", f"{sample}_pseudobulk_counts_barplot.png")
    summary.to_csv(summary_save_path)

    plt.figure(figsize=(10,6))
    sns.barplot(
        data=summary,
        x="cell_type", y="psbulk_counts", hue="CTCFL_expr"
    )

    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Total Raw Counts (psbulk_counts)")
    plt.xlabel("Cell Type")
    plt.title(f"{sample} Pseudobulk raw counts per cell type × CTCFL status")

    plt.tight_layout()
    plt.savefig(summary_barplot_save_path, dpi=300)
    plt.show()

    # --- filter low-quality samples
    dc.pp.filter_samples(pdata, min_cells=5, min_counts=500)
    print("After min_cells/min_counts filter:", pdata.shape)

    # --- QC plot
    if pdata.n_obs > 0:
        dc.pl.obsbar(pdata, y="cell_type", hue="CTCFL_expr", figsize=(6, 3))
        plt.savefig(
            os.path.join(figures_dir, "pseudobulk", f"{sample}_pseudobulk_QC.png"),
            dpi=300, bbox_inches="tight",
        )
        plt.close()

    # --- keep raw counts in .X
    pdata.layers["counts"] = pdata.X.copy()
    dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

    print(f"Remaining pseudobulk samples: {pdata.n_obs}")
    return pdata




def make_deseq_ready(adata, sample, min_gene_total, min_sample_total, n_top_genes):
    """
    Return a new AnnData with integer raw counts in .X for pydeseq2,
    prefiltered to drop weak genes, shallow samples, and (optionally) 
    keep only the top N most variable genes.
    """
    # --- choose counts matrix ---
    if "counts" in adata.layers:
        X = adata.layers["counts"]
    elif adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    # densify if sparse
    if sp.issparse(X):
        X = X.toarray()

    # clean values → non-negative integers
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    X = np.rint(X).clip(min=0).astype(int)

    # --- filter genes ---
    gene_totals = X.sum(axis=0)
    gene_vars = X.var(axis=0)
    keep_genes = (gene_totals >= min_gene_total) & (gene_vars > 0)

    # --- filter samples ---
    sample_totals = X.sum(axis=1)
    keep_samples = sample_totals >= min_sample_total

    # subset once
    X = X[np.ix_(keep_samples, keep_genes)]
    obs = adata.obs.iloc[keep_samples].copy()
    var = adata.var.iloc[keep_genes].copy()

    # --- optionally keep only top N most variable genes ---
    if n_top_genes is not None and X.shape[1] > n_top_genes:
        varscores = X.var(axis=0)
        top_idx = np.argsort(varscores)[-n_top_genes:]
        X = X[:, top_idx]
        var = var.iloc[top_idx]

    print(f"Prefilter: {X.shape[0]} samples × {X.shape[1]} genes remain in sample {sample}"
          f"(≥{min_sample_total} reads per sample, ≥{min_gene_total} counts per gene, var>0).")

    return anndata.AnnData(X=X, obs=obs, var=var, uns=adata.uns.copy())


def pseudobulk_dea(adata_counts, sample, log2_cutoff, pval_cutoff, tables_dir, figures_dir):
    from pydeseq2.dds import DeseqDataSet, DefaultInference
    from pydeseq2.ds import DeseqStats

    os.makedirs(os.path.join(figures_dir, "volcano"), exist_ok=True)
    inference = DefaultInference(n_cpus=8)

    dds = DeseqDataSet(
        adata=adata_counts,
        design_factors=["CTCFL_expr"],
        refit_cooks=True,
        inference=inference,
        n_cpus=12,
    )
    dds.deseq2()

    if set(["CTCFL_pos", "CTCFL_neg"]).issubset(adata_counts.obs["CTCFL_expr"].unique()):
        stat_res = DeseqStats(dds, contrast=["CTCFL_expr", "CTCFL_pos", "CTCFL_neg"], inference=inference)
        stat_res.summary()

        results_df = stat_res.results_df
        results_df.to_csv(os.path.join(tables_dir, "pseudobulk", f"{sample}_deseq2_results.csv"))

        # --- filter significant DEGs ---
        deg = results_df[
            (results_df["padj"] < pval_cutoff) & (results_df["log2FoldChange"].abs() > log2_cutoff)
        ]
        deg.to_csv(os.path.join(tables_dir, "pseudobulk", f"{sample}_deseq2_sig_results.csv"))
        print("Significant DEGs:", deg.shape[0])

        # --- volcano plot ---
        plt.figure(figsize=(7, 6))

        # background
        plt.scatter(results_df["log2FoldChange"], -np.log10(results_df["pvalue"]),
                    c="lightgray", alpha=0.6, s=10, label="All genes")

        # up and down DEGs
        up = deg[deg["log2FoldChange"] > log2_cutoff]
        down = deg[deg["log2FoldChange"] < -log2_cutoff]

        plt.scatter(up["log2FoldChange"], -np.log10(up["pvalue"]),
                    c="red", alpha=0.8, s=12, label="Up DEGs")
        plt.scatter(down["log2FoldChange"], -np.log10(down["pvalue"]),
                    c="blue", alpha=0.8, s=12, label="Down DEGs")

        # thresholds
        plt.axhline(-np.log10(pval_cutoff), color="black", linestyle="--", lw=1)
        plt.axvline(log2_cutoff, color="green", linestyle="--", lw=1)
        plt.axvline(-log2_cutoff, color="green", linestyle="--", lw=1)

        # label top 5 up & down DEGs
        top_hits = pd.concat([
            up.nsmallest(5, "pvalue"),
            down.nsmallest(5, "pvalue")
        ])
        for gene, row in top_hits.iterrows():
            plt.text(row["log2FoldChange"], -np.log10(row["pvalue"]),
                     gene, fontsize=8, ha="right")

        plt.title(f"{sample} Volcano Plot", fontsize=14)
        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10 p-value")
        plt.legend()
        plt.tight_layout()

        plt.savefig(os.path.join(figures_dir, "volcano", f"{sample}_volcano.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

    else:
        print("⚠️ Skipping DESeq2: not enough groups in CTCFL_expr.")


def singlecell_de_volcano(
    adata, sample, figures_dir, tables_dir,
    groupby="CTCFL_expr", method="wilcoxon",
    log2_cutoff=1, pval_cutoff=0.05,
    top_n=10
):
    """Run DE with Scanpy rank_genes_groups and make volcano plot.
    Highlights only the top N up and top N down genes if significant.
    """

    # --- check group sizes
    group_sizes = adata.obs[groupby].value_counts()
    print("Group sizes:", group_sizes.to_dict())
    if group_sizes.min() < 2:
        print(f"⚠️ Skipping {sample}: group {groupby} has fewer than 2 cells.")
        return

    print(f"Running DE for {sample} (groupby={groupby}, method={method})...")

    # --- preprocessing for DE
    if adata.X.max() > 50:   # heuristic: raw counts not normalized
        print("Normalizing and log-transforming data for DE...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # --- DE test
    sc.tl.rank_genes_groups(
        adata, groupby=groupby, method=method,
        use_raw=False, pts=True
    )


    # --- collect results
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    de_dfs = []
    for g in groups:
        df = pd.DataFrame({
            "gene": result["names"][g],
            "pvals": result["pvals"][g],
            "pvals_adj": result["pvals_adj"][g],
            "logfoldchanges": result["logfoldchanges"][g],
        })
        df["group"] = g
        de_dfs.append(df)
    de_df = pd.concat(de_dfs, ignore_index=True)

    # save table
    out_table = os.path.join(tables_dir, "volcano", f"{sample}_{groupby}_DE_results.csv")
    os.makedirs(os.path.dirname(out_table), exist_ok=True)
    de_df.to_csv(out_table, index=False)
    print(f"Saved DE results → {out_table}")

    # --- significance filter
    sig_mask = (de_df["pvals_adj"] < pval_cutoff) & (de_df["logfoldchanges"].abs() > log2_cutoff)
    de_sig = de_df[sig_mask]

    # select top N up and down
    up_genes = (de_sig[de_sig["logfoldchanges"] > 0]
                .nsmallest(top_n, "pvals_adj"))
    down_genes = (de_sig[de_sig["logfoldchanges"] < 0]
                  .nsmallest(top_n, "pvals_adj"))

    # --- volcano plot
    plt.figure(figsize=(7, 6))

    # background (all genes)
    plt.scatter(
        de_df["logfoldchanges"], -np.log10(de_df["pvals"]+1e-300),
        c="lightgray", alpha=0.6, s=10, label="All genes"
    )

    # significant background points (grey)
    plt.scatter(
        de_sig["logfoldchanges"], -np.log10(de_sig["pvals"]+1e-300),
        c="darkgray", alpha=0.7, s=12
    )

    # highlight top up (red) and down (blue)
    plt.scatter(up_genes["logfoldchanges"], -np.log10(up_genes["pvals"]+1e-300),
                c="red", s=25, label=f"Top {top_n} up")
    plt.scatter(down_genes["logfoldchanges"], -np.log10(down_genes["pvals"]+1e-300),
                c="blue", s=25, label=f"Top {top_n} down")

    # label them
    for _, row in pd.concat([up_genes, down_genes]).iterrows():
        plt.text(
            row["logfoldchanges"], -np.log10(row["pvals"]+1e-300),
            row["gene"], fontsize=8, ha="right"
        )

    # thresholds
    plt.axhline(-np.log10(pval_cutoff), color="black", linestyle="--", lw=1)
    plt.axvline(log2_cutoff, color="green", linestyle="--", lw=1)
    plt.axvline(-log2_cutoff, color="green", linestyle="--", lw=1)

    plt.title(f"{sample}: {groupby} DE ({method})", fontsize=14)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 p-value")
    plt.legend()
    plt.tight_layout()

    out_plot = os.path.join(figures_dir, "volcano", f"{sample}_{groupby}_volcano.png")
    os.makedirs(os.path.dirname(out_plot), exist_ok=True)
    plt.savefig(out_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved volcano plot → {out_plot}")
