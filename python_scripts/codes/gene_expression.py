import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np
import scanpy as sc

from config import *

from codes.utils import subset_adata, get_cell_type, add_list_to_config, get_raw_gene_counts, init_color_maps, make_color_map, GENE_COLOR_MAP



init_color_maps(CTCFL_MARKERS + ALL_GENES)
import matplotlib.patches as mpatches
from collections import Counter


def export_all_expression(adata, sample, tables_dir):
    """

    Export the full expression matrix to CSV (genes × cells).
    Priority of source: adata.raw → adata.layers['counts'] → adata.X.
    Also returns the transposed DataFrame (genes as rows).


    """
    if adata.raw is not None:
        ad_use = adata.raw.to_adata()
        X, genes, cells = ad_use.X, ad_use.var_names, ad_use.obs_names
    elif "counts" in adata.layers:
        X, genes, cells = adata.layers["counts"], adata.var_names, adata.obs_names
    else:
        X, genes, cells = adata.X, adata.var_names, adata.obs_names
        if np.issubdtype(X.dtype, np.integer):
            print(f"[{sample}] Using adata.X as raw counts (integer dtype).")
        else:
            print(f"[{sample}] Using adata.X (likely normalized).")

    if sp.issparse(X):
        X = X.toarray()

    df_T = pd.DataFrame(X, index=cells, columns=genes).T

    out_dir = os.path.join(tables_dir, "all_expression")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sample}_all_counts_transpose.csv")
    df_T.to_csv(out_path)
    print(f"Saved all counts → {out_path}")

    return df_T


def export_gene_expression_and_umap(
    adata,
    sample,
    tables_dir,
    figures_dir,
    out_folder,
    important_genes,
    celltype_key="cell_type",
    plot_yes=False,
):
    """
    Export expression tables and generate UMAP plots.

    - Uses only important_genes provided by user.
    - Colors: cell_type_colors and GENE_COLOR_MAP.
    - Binary UMAPs: only top 10 expressed + always CTCFL.
    """

    GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CTCFL_MARKERS + ALL_GENES)

    # --- normalize names once ---
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    all_genes_list = list(adata.var_names)

    # dedupe and normalize important genes
    important_genes = [str(g) for g in (important_genes or [])]
    important_genes = list(dict.fromkeys(important_genes))

    # --- choose data source ---
    if "counts" in adata.layers:
        X, genes, cells = adata.layers["counts"], adata.var_names.astype(str), adata.obs_names.astype(str)
    elif adata.raw is not None:
        ad_use = adata.raw.to_adata()
        X, genes, cells = ad_use.X, ad_use.var_names.astype(str), ad_use.obs_names.astype(str)
    else:
        X, genes, cells = adata.X, adata.var_names.astype(str), adata.obs_names.astype(str)
        print(f"[{sample}] Using adata.X")

    if sp.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, index=cells, columns=genes)

    # --- gene selection ---
    present_genes = [g for g in important_genes if g in all_genes_list]
    missing_genes = [g for g in important_genes if g not in all_genes_list]
    if missing_genes:
        print(f"[{sample}] Skipping missing genes: {missing_genes}")

    if not present_genes:
        print(f"[{sample}] No important genes present, skipping export.")
        return None

    # expression tables
    df_imp = df[present_genes].T
    df_imp["average"] = df_imp.mean(axis=1)
    df_avg = df_imp[["average"]].sort_values(by="average", ascending=False)
    top_genes = df_avg.head(10).index.tolist()
    print("Top genes:", top_genes)

    tables_out = os.path.join(tables_dir, "genes")
    os.makedirs(tables_out, exist_ok=True)
    df_imp.to_csv(os.path.join(tables_out, f"{sample}_genes_cells.csv"))
    df_avg.to_csv(os.path.join(tables_out, f"{sample}_genes.csv"))

    # output dirs
    figs_out = os.path.join(figures_dir, out_folder)
    summary_out = os.path.join(figs_out, "summary_marked_genes")
    binary_out = os.path.join(figs_out, "binary_marked_genes")
    expr_out = os.path.join(figs_out, "expression_marked_genes")
    os.makedirs(summary_out, exist_ok=True)
    os.makedirs(binary_out, exist_ok=True)
    os.makedirs(expr_out, exist_ok=True)

    # --- ensure UMAP embedding exists ---
    if "X_umap" not in adata.obsm:
        print(f"[{sample}] No UMAP embedding found → computing one.")
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    # 1) Cell type UMAP
    cats = adata.obs[celltype_key].astype(str).unique().tolist()

    # build a palette that only includes present categories
    palette = {c: CELL_TYPE_COLORS.get(c, "#d3d3d3") for c in cats}

    fig = sc.pl.umap(
        adata,
        color=celltype_key,
        frameon=False,
        show=False,
        return_fig=True,
        title=f"{sample}\nPredicted Cell Type",
        palette=palette,  # safe palette
    )

    

    fig.set_size_inches(6, 6)
    fig.savefig(os.path.join(summary_out, f"celltype_{sample}_umap.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 2) Combined categorical UMAP (top 5 genes, skip None in legend)
    if present_genes:
        mat = adata[:, present_genes].layers["counts"] if "counts" in adata.layers else adata[:, present_genes].X
        mat = mat.toarray() if sp.issparse(mat) else mat
        max_idx = np.argmax(mat, axis=1)
        labels = [present_genes[i] if mat[j, i] > 0 else "None"
                  for j, i in enumerate(max_idx)]
        adata.obs["CTCFL_marker"] = labels

        counts = Counter(labels)
        top5 = [g for g, _ in counts.most_common() if g != "None"][:5]

        color_map = {g: GENE_COLOR_MAP.get(g, "#d3d3d3") for g in present_genes}
        color_map["None"] = "#d3d3d3"
        for g in color_map:
            if g not in top5 and g != "None":
                color_map[g] = "#d3d3d3"

        fig = sc.pl.umap(
            adata,
            color="CTCFL_marker",
            palette=color_map,
            frameon=False,
            show=False,
            return_fig=True,
            title=f"{sample} — Overall Top Expressed Gene/Cell \n (Top 5 Per Sample)"
        )

        legend_handles = [
            mpatches.Patch(color=color_map[g], label=f"{g} ({counts[g]})")
            for g in top5
        ]
        for ax in fig.axes:
            ax.legend(handles=legend_handles, loc="best", frameon=False)

        fig.set_size_inches(7, 7)
        fig.savefig(os.path.join(summary_out, f"{sample}_top5_genes_umap.png"),
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 3) Binary UMAPs: only top 10 + always CTCFL
    if plot_yes:
        binary_genes = list(top_genes)

        # ensure CTCFL is always included
        ctcfl_matches = [g for g in all_genes_list if g.upper() == "CTCFL"]
        if ctcfl_matches:
            ctcfl_gene = ctcfl_matches[0]
            if ctcfl_gene not in binary_genes:
                binary_genes.append(ctcfl_gene)
            print(f"✅ Ensured CTCFL included as {ctcfl_gene}")
        else:
            print("⚠️ CTCFL not found in dataset, adding placeholder all-neg map")
            ctcfl_gene = "CTCFL"
            binary_genes.append(ctcfl_gene)

        for g in binary_genes:
            if g not in df.columns:
                print(f"[{sample}] Skipping {g} (not in dataset)")
                continue

            expr = df[g].values
            pos_mask = expr > 0
            n_pos, n_total = int(pos_mask.sum()), expr.shape[0]

            adata.obs[f"{g}_binary"] = pd.Categorical(
                np.where(pos_mask, "pos", "neg"),
                categories=["neg", "pos"]
            )

            if g not in GENE_COLOR_MAP:
                make_color_map([g])  # assign a new color

            pos_color = GENE_COLOR_MAP[g]
            palette = {"neg": "#d3d3d3", "pos": pos_color if n_pos > 0 else "#d3d3d3"}

            print(f"[{sample}] {g}: {n_pos}/{n_total} positive cells")

            fig = sc.pl.umap(
                adata,
                color=f"{g}_binary",
                frameon=False,
                palette=palette,
                show=False,
                return_fig=True,
                title=f"{sample}\n{g} ({n_pos}/{n_total} cells positive)"
            )
            fig.set_size_inches(6, 6)
            fig.savefig(os.path.join(binary_out, f"{sample}_{g}_binary_umap.png"),
                        dpi=300, bbox_inches="tight")
            plt.close(fig)

            # Expression UMAP
            expr_proc = df[g].values
            expr_masked = np.where(expr_proc > 0, expr_proc, np.nan)
            adata.obs[f"{g}_expr_umap"] = expr_masked

            fig = sc.pl.umap(
                adata,
                color=f"{g}_expr_umap",
                cmap="RdYlGn",
                na_color="#d3d3d3",
                frameon=False,
                show=False,
                return_fig=True,
                title=f"{sample}\n{g} Expression"
            )
            fig.set_size_inches(6, 6)
            fig.savefig(os.path.join(binary_out, f"{sample}_{g}_expr_umap.png"),
                        dpi=300, bbox_inches="tight")
            plt.close(fig)

    return df_avg















def export_gene_split(adata, sample, gene, tables_dir):
    """
Split cells into gene+ (expr > 0) and gene– (expr = 0).
Saves:
- Single-cell expression matrices for positive and negative cells.
- Cell type–aggregated average expression for positive and negative cells.
"""


    # --- pick expression source ---
    if adata.raw is not None and gene in adata.raw.var_names:
        Xg = adata.raw[:, gene].X
        X_all = adata.raw.X
        genes = adata.raw.var_names
        print(f"[{sample}] Using adata.raw for {gene}")
    elif "counts" in adata.layers and gene in adata.var_names:
        Xg = adata[:, gene].layers["counts"]
        X_all = adata.layers["counts"]
        genes = adata.var_names
        print(f"[{sample}] Using adata.layers['counts'] for {gene}")
    elif gene in adata.var_names:
        Xg = adata[:, gene].X
        X_all = adata.X
        genes = adata.var_names
        print(f"[{sample}] Using adata.X for {gene}")
    else:
        raise KeyError(f"Gene '{gene}' not found in adata.var_names or adata.raw.var_names.")

    # --- expression vector ---
    gene_expr = Xg.toarray().ravel() if sp.issparse(Xg) else np.ravel(Xg)
    mask_gene = gene_expr > 0
    adata.obs[f"{gene}_expr"] = np.where(mask_gene, f"{gene}_pos", f"{gene}_neg")
    adata.obs[f"{gene}_expr"] = adata.obs[f"{gene}_expr"].astype("category")

    print(f"\n{mask_gene.sum()} / {adata.n_obs} cells are {gene}_pos")

    # --- full matrix ---
    if sp.issparse(X_all):
        X_all = X_all.toarray()
    df = pd.DataFrame(X_all, index=adata.obs_names, columns=genes)

    # --- split single-cell ---
    pos_cells = df.index[mask_gene]
    neg_cells = df.index[~mask_gene]
    df_pos_T, df_neg_T = df.loc[pos_cells].T, df.loc[neg_cells].T

    # --- save ---
    out_dir = os.path.join(tables_dir, f"{gene}_split")
    os.makedirs(out_dir, exist_ok=True)

    if not df_pos_T.empty:
        out_path = os.path.join(out_dir, f"{sample}_{gene}_pos_counts_singlecell.csv")
        df_pos_T.to_csv(out_path)
        print(f"Saved {gene}+ single-cell counts → {out_path}")
    else:
        print(f"[{sample}] No {gene}+ cells found.")

    if not df_neg_T.empty:
        out_path = os.path.join(out_dir, f"{sample}_{gene}_neg_counts_singlecell.csv")
        df_neg_T.to_csv(out_path)
        print(f"Saved {gene}- single-cell counts → {out_path}")
    else:
        print(f"[{sample}] No {gene}- cells found.")

    # --- aggregated by cell_type ---
    if "cell_type" in adata.obs.columns:
        df["cell_type"] = adata.obs["cell_type"]

        agg_pos = df.loc[pos_cells].groupby("cell_type").mean().T
        agg_neg = df.loc[neg_cells].groupby("cell_type").mean().T

        if not agg_pos.empty:
            out_path = os.path.join(out_dir, f"{sample}_{gene}_pos_counts_by_celltype.csv")
            agg_pos.to_csv(out_path)
            print(f"Saved {gene}+ aggregated counts → {out_path}")
        if not agg_neg.empty:
            out_path = os.path.join(out_dir, f"{sample}_{gene}_neg_counts_by_celltype.csv")
            agg_neg.to_csv(out_path)
            print(f"Saved {gene}- aggregated counts → {out_path}")
    else:
        print(f"[{sample}] No 'cell_type' column found → skipping aggregation.")

    return df_pos_T, df_neg_T




    

def plot_gene_expression(df, df_T, sample, gene, figures_dir, top_n=20):
    """
Plot gene expression summary for a given marker gene.

1. Heatmap of top-N most highly expressed genes in gene+ cells (single-cell level).
2. Dotplot of mean expression by cell type (aggregated).
Saves both plots as PNGs in figures_dir/gene/.
"""

    out_dir = os.path.join(figures_dir, gene.lower())
    os.makedirs(out_dir, exist_ok=True)

    # --- Heatmap ---
    top_genes = df_T.mean(axis=1).sort_values(ascending=False).head(top_n).index
    mat = df_T.loc[top_genes]

    # Clean col labels (remove sample + gene pos/neg suffixes)
    clean_cols = (
        mat.columns
        .str.replace(f"{sample}_", "", regex=False)
        .str.replace(f"_{gene}_pos", "", regex=False)
        .str.replace(f"_{gene}_neg", "", regex=False)
    )
    mat.columns = clean_cols

    plt.figure(figsize=(10, max(6, top_n * 0.4)))
    sns.heatmap(mat, cmap="viridis", cbar_kws={"label": "Raw counts"})
    plt.title(f"{sample}: Top {top_n} {gene}+ genes (single-cell)", fontsize=14)
    plt.xlabel("Cells"); plt.ylabel("Genes")
    plt.xticks(rotation=90)
    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{sample}_{gene}_heatmap.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close()
    print(f"Saved heatmap → {out_path}")

    # --- Dotplot ---
    if "cell_type" not in df.columns:
        print(f"[{sample}] No 'cell_type' in df, skipping dotplot.")
        return

    expr, meta = df.drop(columns=["cell_type"]), df["cell_type"]

    # Clean up meta labels too
    meta_clean = (
        meta.astype(str)
        .str.replace(f"{sample}_", "", regex=False)
        .str.replace(f"_{gene}_pos", "", regex=False)
        .str.replace(f"_{gene}_neg", "", regex=False)
    )

    results = []
    for ct in meta_clean.unique():
        mean_expr = expr[meta_clean == ct].mean(axis=0)
        results.append(pd.DataFrame({
            "cell_type": ct, "gene": mean_expr.index, "mean_expr": mean_expr.values
        }))
    if not results:
        print(f"[{sample}] No groups for dotplot.")
        return
    plot_df = pd.concat(results)

    top_genes = (plot_df.groupby("gene")["mean_expr"]
                 .mean().sort_values(ascending=False).head(top_n).index)
    plot_df = plot_df[plot_df["gene"].isin(top_genes)]

    plt.figure(figsize=(8, max(4, top_n * 0.3)))
    sns.scatterplot(data=plot_df, x="cell_type", y="gene",
                    hue="mean_expr", palette="viridis",
                    s=80, edgecolor="black", alpha=0.8)
    plt.title(f"{sample}: Top {top_n} {gene}+ genes (by cell type)", fontsize=14)
    plt.xlabel("Cell type"); plt.ylabel("Genes")
    plt.xticks(rotation=90); sns.despine(); plt.grid(False)
    plt.legend(title="Mean expr", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{sample}_{gene}_dotplot.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close()
    print(f"Saved dotplot → {out_path}")

def find_coexpressed_genes(
    adata,
    sample,
    target_gene,
    tables_dir,
    figures_dir,
    out_folder,
):
    """
    Find genes most co-expressed with a given target gene in the
    dataset (prioritize counts).

    - Computes correlation of all genes vs target_gene.
    - Saves top correlated genes (raw counts) to CSV.
    - Plots UMAPs:
        1. Categorical (cells labeled by top co-expressed gene)
        2. Binary per top_n_umap genes
    """

    # --- normalize names once ---
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)

    # universal colors
    GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CTCFL_MARKERS + ALL_GENES)

    # --- choose data source ---
    if adata.raw is not None:
        ad_use = adata.raw.to_adata()
        X, genes, cells = ad_use.X, ad_use.var_names.astype(str), ad_use.obs_names.astype(str)
    elif "counts" in adata.layers:
        X, genes, cells = adata.layers["counts"], adata.var_names.astype(str), adata.obs_names.astype(str)
    else:
        X, genes, cells = adata.X, adata.var_names.astype(str), adata.obs_names.astype(str)
        print(f"[{sample}] Using adata.X")

    if sp.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, index=cells, columns=genes)

    # --- target gene check (case-insensitive) ---
    matches = [g for g in df.columns if g.upper() == target_gene.upper()]
    if not matches:
        print(f"[{sample}] Target gene {target_gene} not found in adata.var_names")
        return None

    gene = matches[0]
    gene_expr = df[gene].values
    binary_expr = (gene_expr > 0.1).astype(int)
    n_pos = int(binary_expr.sum())

    if n_pos == 0:
        print(f"[{sample}] Skipping: {target_gene} not expressed in any cell.")
        return None
    else:
        print(f"[{sample}] {target_gene}: {n_pos}/{len(binary_expr)} positive cells")

    # --- processed correlations ---
    X_proc = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    df_proc = pd.DataFrame(X_proc, index=adata.obs_names, columns=adata.var_names)

    matches_proc = [g for g in df_proc.columns if g.upper() == target_gene.upper()]
    if not matches_proc:
        print(f"[{sample}] Target gene {target_gene} not found in processed var_names")
        return None

    target_gene_proc = matches_proc[0]
    target_expr = df_proc[target_gene_proc].values

    corrs = df_proc.corrwith(pd.Series(target_expr, index=df_proc.index))
    corrs = corrs.drop(target_gene_proc).sort_values(ascending=False)

    # --- thresholds ---
    top_csv_genes = corrs[corrs > 0.2].index.tolist()
    pos_counts = (df_proc[top_csv_genes] > 0.1).sum()
    top_umap_genes = pos_counts[pos_counts > 4].index.tolist()

    print(f"[{sample}] Genes with corr > 0.2 vs {target_gene}: {top_csv_genes[:10]} ...")
    print(f"[{sample}] Genes with >5 positive counts for UMAPs: {top_umap_genes[:10]}")

    if not top_umap_genes:
        print(f"[{sample}] No genes passed thresholds → skipping UMAPs and moving on.")
        return None

    print(f"[{sample}] Top Correlated vs {target_gene}: {top_csv_genes[:10]} ...")

    # ----- output dirs -----
    tables_out = os.path.join(tables_dir, out_folder, "coexpressed")
    figs_out = os.path.join(figures_dir, out_folder, "coexpressed")
    os.makedirs(tables_out, exist_ok=True)
    os.makedirs(figs_out, exist_ok=True)

    # save correlations
    corrs.head(len(top_csv_genes)).to_csv(
        os.path.join(tables_out, f"{sample}_{target_gene}_top_correlations_processed_data.csv")
    )

    # save raw counts
    valid_top_csv_genes = [g for g in top_csv_genes if g in df.columns]
    if not valid_top_csv_genes:
        print(f"[{sample}] No valid top genes found in raw dataset, skipping CSV export.")
        return adata

    df_sel_csv = df[valid_top_csv_genes]
    df_sel_csv.T.to_csv(os.path.join(tables_out, f"{sample}_{target_gene}_top_raw_counts.csv"))

    # save processed counts
    valid_top_proc_genes = [g for g in top_csv_genes if g in df_proc.columns]
    if valid_top_proc_genes:
        df_norm_csv = df_proc[valid_top_proc_genes]
        df_norm_csv.T.to_csv(os.path.join(tables_out, f"{sample}_{target_gene}_top_processed_data.csv"))

    # --- BORIS+ restriction ---
    if "CTCFL_expr" in adata.obs:
        boris_mask = adata.obs["CTCFL_expr"] == "CTCFL_pos"
        df_boris = df_proc.loc[boris_mask, valid_top_proc_genes]

        if not df_boris.empty:
            out_path = os.path.join(
                tables_out, f"{sample}_{target_gene}_top_processed_BORISpos.csv"
            )
            df_boris.T.to_csv(out_path)
            print(f"[{sample}] Saved BORIS+ cells × top correlated genes → {out_path}")
        else:
            print(f"[{sample}] No BORIS+ cells found, skipping BORIS+ table")

    # ----- ensure UMAP embedding exists -----
    if "X_umap" not in adata.obsm:
        print(f"[{sample}] No UMAP embedding found → computing one.")
        sc.pp.neighbors(adata)

        if adata.n_obs == 0:
            print(f"[{sample}] Empty AnnData, skipping.")
            return None

        sc.tl.umap(adata)

    # ----- categorical UMAP -----
    valid_umap_genes = [g for g in top_umap_genes if g in df.columns]
    if not valid_umap_genes:
        print(f"[{sample}] No valid UMAP genes found in raw dataset → skipping categorical UMAP.")
        return adata

    mat = df[valid_umap_genes].values
    max_idx = np.argmax(mat, axis=1)
    labels = [valid_umap_genes[i] if mat[j, i] > 0 else "None" for j, i in enumerate(max_idx)]
    adata.obs[f"{target_gene}_coexpr_marker"] = labels

    color_map = {g: GENE_COLOR_MAP.get(g, "#d3d3d3") for g in valid_umap_genes}
    color_map["None"] = "#d3d3d3"

    fig = sc.pl.umap(
        adata,
        color=f"{target_gene}_coexpr_marker",
        palette=color_map,
        frameon=False,
        show=False,
        return_fig=True,
        title=f"{sample} — Top co-expressed with {target_gene}"
    )
    fig.set_size_inches(7, 7)
    fig.savefig(
        os.path.join(figs_out, f"{sample}_{target_gene}_categorical_umap.png"),
        dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # ----- binary UMAPs -----
    valid_umap_genes.append(target_gene)
    for g in valid_umap_genes:
        if g not in df.columns:
            continue
        expr = df[g].values
        binary_expr = (expr >= 1).astype(int)

        n_pos = int(binary_expr.sum())
        n_total = len(binary_expr)

        adata.obs[f"{g}_binary_coexpr"] = pd.Categorical(
            np.where(binary_expr.astype(bool), "pos", "neg"),
            categories=["neg", "pos"]
        )

        if g not in GENE_COLOR_MAP:
            make_color_map([g])  # assign new color

        pos_color = GENE_COLOR_MAP[g]
        palette = {"neg": "#d3d3d3", "pos": pos_color}

        fig = sc.pl.umap(
            adata,
            color=f"{g}_binary_coexpr",
            palette=palette,
            frameon=False,
            show=False,
            return_fig=True,
            title=f"{sample}\n{g} ({n_pos}/{n_total} cells positive)"
        )
        fig.set_size_inches(7, 7)
        fig.savefig(
            os.path.join(figs_out, f"{sample}_{g}_coexpressed_binary_umap.png"),
            dpi=300, bbox_inches="tight"
        )
        plt.close(fig)

    return adata
