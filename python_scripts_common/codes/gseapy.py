
# projects/ovarian_cancer/python_scripts/codes/gseapy.py

#### Libraries ###
import os
import glob
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import gseapy as gp

#### Custom ####
from config import *

### Functions ###

def load_all_correlations(
    tables_dir=TABLES_DIR,
    gene_filter_regex=r"^(AC|AP|LINC)",
    drop_gene="CTCFL",
):
    corr_dir = os.path.join(tables_dir, "DE_analysis")
    files = glob.glob(os.path.join(corr_dir, "*_CTCFL_correlations.csv"))

    dfs = []
    for f in files:
        sample = os.path.basename(f).split("_")[0]
        df = pd.read_csv(f)
        df["sample"] = sample
        dfs.append(df)

    all_corr = pd.concat(dfs, ignore_index=True)

    all_corr = all_corr[
        (all_corr["gene"] != drop_gene)
        & (~all_corr["gene"].str.match(gene_filter_regex))
    ]

    return all_corr


def split_pos_neg_programs(
    all_corr,
    corr_thresh=0.3,
    min_samples=3,
):
    corr_mat = all_corr.pivot_table(
        index="gene",
        columns="sample",
        values="correlation"
    )

    mean_corr = corr_mat.mean(axis=1)
    n_samples = corr_mat.notna().sum(axis=1)

    pos_genes = corr_mat.loc[
        (mean_corr > corr_thresh) & (n_samples >= min_samples)
    ]

    neg_genes = corr_mat.loc[
        (mean_corr < -corr_thresh) & (n_samples >= min_samples)
    ]

    pos_genes = pos_genes.assign(mean_corr=mean_corr)
    pos_genes = pos_genes.sort_values("mean_corr", ascending=False)

    neg_genes = neg_genes.assign(mean_corr=mean_corr)
    neg_genes = neg_genes.sort_values("mean_corr")

    return pos_genes, neg_genes


def run_enrichr(
    gene_list,
    gene_sets=None,
):
    if gene_sets is None:
        gene_sets = [
            "Reactome_2022",
            "GO_Biological_Process_2023",
            "MSigDB_Hallmark_2020",
        ]

    enr = gp.enrichr(
        gene_list=gene_list,
        organism="Human",
        gene_sets=gene_sets,
        outdir=None,
    )

    return enr.results


def plot_enrichment_bar(
    enr_df,
    title,
    out_path,
    top_n=10,
    color="#4C72B0",
):
    df = (
        enr_df
        .sort_values("Adjusted P-value")
        .head(top_n)
        .iloc[::-1]
        .copy()
    )

    df["score"] = -np.log10(df["Adjusted P-value"])

    plt.figure(figsize=(6, 4))
    plt.barh(df["Term"], df["score"], color=color)
    plt.xlabel("-log10(adj p-value)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_pos_neg_comparison(
    enr_pos,
    enr_neg,
    out_path,
    top_n=8,
):
    def prep(df, label):
        d = (
            df
            .sort_values("Adjusted P-value")
            .head(top_n)
            .copy()
        )
        d["group"] = label
        d["score"] = -np.log10(d["Adjusted P-value"])
        return d[["Term", "score", "group"]]

    combo = pd.concat([
        prep(enr_pos, "Positive"),
        prep(enr_neg, "Negative"),
    ])

    plt.figure(figsize=(7, 5))
    sns.barplot(
        data=combo,
        x="score",
        y="Term",
        hue="group",
        palette={
            "Positive": "#D55E00",
            "Negative": "#0072B2",
        }
    )
    plt.xlabel("-log10(adj p-value)")
    plt.ylabel("")
    plt.title("CTCFL pathway programs")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
