#### Libraries ####
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import scipy.sparse as sp
from collections import Counter

# Custom
from setup import init_color_maps, make_color_map, GENE_COLOR_MAP, CELL_TYPE_COLORS

from config import *


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
    - Colors: CELL_TYPE_COLORS and GENE_COLOR_MAP.
    - Binary UMAPs: top 10 expressed genes + always CTCFL.
    """

    # --- init palettes ---
    GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CTCFL_MARKERS + ALL_GENES)

    # normalize names
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    all_genes = list(adata.var_names)

    # dedupe important genes
    important_genes = list(dict.fromkeys(map(str, important_genes or [])))

    # --- choose data source ---
    if "counts" in adata.layers:
        X, genes, cells = adata.layers["counts"], adata.var_names, adata.obs_names
    elif adata.raw is not None:
        ad_use = adata.raw.to_adata()
        X, genes, cells = ad_use.X, ad_use.var_names, ad_use.obs_names
    else:
        print(f"[{sample}] Using adata.X")
        X, genes, cells = adata.X, adata.var_names, adata.obs_names

    if sp.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, index=cells, columns=genes)

    # --- gene selection ---
    present = [g for g in important_genes if g in all_genes]
    missing = [g for g in important_genes if g not in all_genes]
    if missing:
        print(f"[{sample}] Skipping missing genes: {missing}")
    if not present:
        print(f"[{sample}] No important genes present, skipping export.")
        return None

    # expression stats
    df_imp = df[present].T
    df_imp["average"] = df_imp.mean(axis=1)
    df_avg = df_imp[["average"]].sort_values(by="average", ascending=False)
    top_genes = df_avg.head(10).index.tolist()
    print("Top genes:", top_genes)

    # --- output dirs ---
    tables_out = os.path.join(tables_dir, "genes")
    figs_out = os.path.join(figures_dir, out_folder)
    summary_out = os.path.join(figs_out, "summary_marked_genes")
    marked_out = os.path.join(figs_out, "marked_genes")  # merged folder
    for d in [tables_out, summary_out, marked_out]:
        os.makedirs(d, exist_ok=True)

    # save tables
    df_imp.to_csv(os.path.join(tables_out, f"{sample}_genes_cells.csv"))
    df_avg.to_csv(os.path.join(tables_out, f"{sample}_genes.csv"))

    # --- ensure UMAP ---
    if "X_umap" not in adata.obsm:
        print(f"[{sample}] No UMAP embedding found → computing one.")
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    # 1) Cell type UMAP
    cats = adata.obs[celltype_key].astype(str).unique().tolist()
    palette = {c: CELL_TYPE_COLORS.get(c, "#d3d3d3") for c in cats}
    fig = sc.pl.umap(
        adata,
        color=celltype_key,
        frameon=False,
        show=False,
        return_fig=True,
        title=f"{sample}\nPredicted Cell Type",
        palette=palette,
    )
    fig.set_size_inches(6, 6)
    fig.savefig(os.path.join(summary_out, f"celltype_{sample}_umap.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)



    # 2) Top-5 combined categorical UMAP
    mat = adata[:, present].layers["counts"] if "counts" in adata.layers else adata[:, present].X
    mat = mat.toarray() if sp.issparse(mat) else mat
    max_idx = np.argmax(mat, axis=1)
    labels = [present[i] if mat[j, i] > 0 else "None" for j, i in enumerate(max_idx)]
    adata.obs["CTCFL_marker"] = labels

    counts = Counter(labels)
    top5 = [g for g, _ in counts.most_common() if g != "None"][:5]

    # start with defaults
    color_map = {}
    for g in present:
        if g in top5:
            color_map[g] = GENE_COLOR_MAP.get(g, "#d3d3d3")
        else:
            color_map[g] = "#d3d3d3"

    # always add "None"
    color_map["None"] = "#d3d3d3"

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

    # 3) Binary + expression UMAPs
    if plot_yes:
        binary_genes = top_genes[:]
        if "CTCFL" not in [g.upper() for g in binary_genes]:
            ctcfl = next((g for g in all_genes if g.upper() == "CTCFL"), None)
            if ctcfl:
                binary_genes.append(ctcfl)
                print(f"✅ Ensured CTCFL included as {ctcfl}")
            else:
                print("⚠️ CTCFL not found, adding placeholder map")
                binary_genes.append("CTCFL")

        for g in binary_genes:
            if g not in df.columns:
                print(f"[{sample}] Skipping {g} (not in dataset)")
                continue

            expr = df[g].values
            pos_mask = expr > 0
            n_pos, n_total = int(pos_mask.sum()), expr.shape[0]

            adata.obs[f"{g}_binary"] = pd.Categorical(
                np.where(pos_mask, "pos", "neg"), categories=["neg", "pos"]
            )
            if g not in GENE_COLOR_MAP:
                make_color_map([g])

            palette = {"neg": "#d3d3d3", "pos": GENE_COLOR_MAP[g] if n_pos > 0 else "#d3d3d3"}
            print(f"[{sample}] {g}: {n_pos}/{n_total} positive cells")

            # binary map
            fig = sc.pl.umap(
                adata,
                color=f"{g}_binary",
                frameon=False,
                palette=palette,
                show=False,
                return_fig=True,
                title=f"{sample}\n{g} ({n_pos}/{n_total} cells pos)"
            )
            fig.set_size_inches(6, 6)
            fig.savefig(os.path.join(marked_out, f"{sample}_{g}_binary.png"), dpi=300)
            plt.close(fig)

            # expression map
            adata.obs[f"{g}_expr_umap"] = np.where(expr > 0, expr, np.nan)
            fig = sc.pl.umap(
                adata,
                color=f"{g}_expr_umap",
                cmap="RdYlGn",
                na_color="#d3d3d3",
                frameon=False,
                show=False,
                return_fig=True,
                title=f"{sample}\n{g} ({n_pos}/{n_total} cells pos)"
            )
            fig.set_size_inches(6, 6)
            fig.savefig(os.path.join(marked_out, f"{sample}_{g}_expr.png"), dpi=300)
            plt.close(fig)

    return df_avg











