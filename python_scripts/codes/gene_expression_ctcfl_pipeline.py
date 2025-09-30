import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np
import scanpy as sc

from config import *


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
    ctcfl_colors=None, cell_type_colors=None  #
):
    """
    Export expression tables and generate UMAP plots for important genes.

    - Saves gene × cell expression tables and average expression summaries.
    - Plots:
        1. UMAP colored by cell type.
        2. UMAP showing top 5 most co-expressed marker genes.
        3. (Optional) binary UMAPs per gene (positive vs negative cells).
"""

    # ----- choose data source -----
    if adata.raw is not None:
        ad_use = adata.raw.to_adata()
        X, genes, cells = ad_use.X, ad_use.var_names, ad_use.obs_names
    elif "counts" in adata.layers:
        X, genes, cells = adata.layers["counts"], adata.var_names, adata.obs_names
    else:
        X, genes, cells = adata.X, adata.var_names, adata.obs_names
        print(f"[{sample}] Using adata.X")

    if sp.issparse(X):
        X = X.toarray()

    # ----- build dataframe -----
    df = pd.DataFrame(X, index=cells, columns=genes)

    present_genes = [g for g in important_genes if g in adata.var_names]
    missing_genes = [g for g in important_genes if g not in adata.var_names]

    if missing_genes:
        print(f"[{sample}] Skipping missing genes: {missing_genes}")
    df_imp = df[present_genes].T

    df_imp["average"] = df_imp.mean(axis=1)
    df_avg = df_imp[["average"]].sort_values(by="average", ascending=False)
    top_genes = df_avg.head(10).index.tolist()
    print("Top important genes:", top_genes)

    # ----- Save tables -----
    tables_out = os.path.join(tables_dir, "genes")
    os.makedirs(tables_out, exist_ok=True)
    table_path = os.path.join(tables_out, f"{sample}_important_genes_cells.csv")
    table_path_avg = os.path.join(tables_out, f"{sample}_important_genes.csv")
    df_imp.to_csv(table_path)
    df_avg.to_csv(table_path_avg)
    print(f"Saved important genes ({len(present_genes)}) → {table_path}")

    # ----- Output folders -----
    figs_out = os.path.join(figures_dir, out_folder)
    os.makedirs(figs_out, exist_ok=True)
    summary_out = os.path.join(figs_out, "summary_marked_genes")
    binary_out = os.path.join(figs_out, "binary_marked_genes")
    os.makedirs(binary_out, exist_ok=True)
    os.makedirs(summary_out, exist_ok=True)

    # 1) Cell type UMAP
    fig = sc.pl.umap(
        adata,
        color=celltype_key,
        frameon=False,
        show=False,
        return_fig=True, title=f"{sample}\n Predicted Cell Type",
        palette=cell_type_colors   #
    )
    fig.set_size_inches(6, 6)
    fig_path = os.path.join(summary_out, f"celltype_{sample}_umap.png")
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {fig_path}")

    # 3) Combined categorical UMAP
        # 3) Combined categorical UMAP
    if len(present_genes) > 0:
        mat = adata[:, present_genes].X
        mat = mat.toarray() if sp.issparse(mat) else mat

        import numpy as np, matplotlib.patches as mpatches
        from collections import Counter

        max_idx = np.argmax(mat, axis=1)
        labels = [present_genes[i] if mat[j, i] > 0 else "None"
                  for j, i in enumerate(max_idx)]
        adata.obs["CTCFL_marker"] = labels

        counts = Counter(labels)
        top_genes = [g for g, _ in counts.most_common() if g not in ["None", "CTCFL"]][:5]

        color_map = {g: ctcfl_colors.get(g, "#d3d3d3") for g in present_genes}
        color_map["None"] = "#d3d3d3"
        if "CTCFL" in color_map:
            color_map["CTCFL"] = "#d3d3d3"

        for g in color_map:
            if g not in top_genes and g not in ["None", "CTCFL"]:
                color_map[g] = "#d3d3d3"

        fig = sc.pl.umap(
            adata,
            color="CTCFL_marker",
            palette=color_map,
            frameon=False,
            show=False,
            return_fig=True,
            title=f"{sample} — top 5 co-expressed markers per cell"
        )

        legend_handles = [
            mpatches.Patch(color=ctcfl_colors.get(g, "#d3d3d3"), label=f"{g}")
            for g in top_genes
        ]
        for ax in fig.axes:
            ax.legend(handles=legend_handles, loc="best", frameon=False)

        fig.set_size_inches(7, 7)
        fig_path = os.path.join(summary_out, f"{sample}_top5_markers_umap.png")
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {fig_path}")

    # 2) Binary UMAPs (per gene) — restored ✅
    if plot_yes:
        import matplotlib.patches as mpatches
        for g in present_genes:
            if "counts" in adata.layers:
                mat, genes = adata.layers["counts"], adata.var_names
            elif adata.raw is not None:
                mat, genes = adata.raw.X, adata.raw.var_names
            else:
                mat, genes = adata.X, adata.var_names

            if g not in genes:
                print(f"⚠️ {g} not found, skipping")
                continue

            idx = list(genes).index(g)
            expr = mat[:, idx]
            expr = expr.toarray().ravel() if sp.issparse(expr) else np.ravel(expr)

            pos_mask = expr > 0
            n_pos, n_total = int(pos_mask.sum()), expr.shape[0]

            if n_pos == 0:
                adata.obs[f"{g}_binary"] = pd.Categorical(["neg"] * n_total, categories=["neg"])
                palette = {"neg": "#d3d3d3"}
                legend_handles = [mpatches.Patch(color="#d3d3d3", label=f"{g}_neg (0/{n_total})")]
            else:
                adata.obs[f"{g}_binary"] = pd.Categorical(
                    np.where(pos_mask, "pos", "neg"), categories=["neg", "pos"]
                )
                # ✅ Always use predefined BORIS/marker color
                pos_color = ctcfl_colors.get(g, "#ff7f0e")  
                palette = {"neg": "#d3d3d3", "pos": pos_color}
                legend_handles = [
                    mpatches.Patch(color="#d3d3d3", label=f"{g}_neg"),
                    mpatches.Patch(color=pos_color, label=f"{g}_pos ({n_pos}/{n_total})"),
                ]

            fig = sc.pl.umap(
                adata,
                color=f"{g}_binary",
                frameon=False,
                palette=palette,
                show=False,
                return_fig=True,
                title=f"{sample}\n{g} ({n_pos}/{n_total} cells positive)"
            )

            for ax in fig.axes:
                if hasattr(ax, "collections"):
                    for coll in ax.collections:
                        if hasattr(coll, "colorbar") and coll.colorbar:
                            coll.colorbar.remove()
                ax.legend(handles=legend_handles, loc="best", frameon=False)

            fig.set_size_inches(6, 6)
            fig_path = os.path.join(binary_out, f"{sample}_{g}_binary_umap.png")
            fig.savefig(fig_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"Saved → {fig_path}")




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


def find_ctcfl_signature(adata, sample, tables_dir, gene = "CTCFL"):
    if os.path.exists(os.path.join(tables_dir, "gene", f"{sample}_{gene}_neg.csv")):
        df_pos_T = os.path.join(tables_dir, "gene", f"{sample}_{gene}_pos.csv")
        df_neg_T = os.path.join(tables_dir, "gene", f"{sample}_{gene}_neg.csv")
    else:
        df_pos_T, df_neg_T = export_gene_split(adata, sample, gene, tables_dir)
    
    print(df_pos_T.head())

adata = sc.read_h5ad("working_adata_tumors/T76.h5ad")
sample = "T76"
tables_dir = "results/tables"

find_ctcfl_signature(adata, sample, tables_dir, gene = "CTCFL")



