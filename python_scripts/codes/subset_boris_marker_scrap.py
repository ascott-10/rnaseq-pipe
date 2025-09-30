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
    Same as before, but restricted to boris_marker cells.
    """

    # --- subset boris_marker cells only ---
    if "BORIS marker" not in adata.obs[celltype_key].unique():
        print(f"[{sample}] No boris_marker cells found, skipping.")
        return None
    adata = adata[adata.obs[celltype_key] == "boris_marker"].copy()

    GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CTCFL_MARKERS + ALL_GENES)

    # normalize names
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    all_genes_list = list(adata.var_names)

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
    if not present_genes:
        print(f"[{sample}] No important genes present, skipping export.")
        return None

    df_imp = df[present_genes].T
    df_imp["average"] = df_imp.mean(axis=1)
    df_avg = df_imp[["average"]].sort_values(by="average", ascending=False)
    top_genes = df_avg.head(10).index.tolist()

    # save tables
    tables_out = os.path.join(tables_dir, "genes")
    os.makedirs(tables_out, exist_ok=True)
    df_imp.to_csv(os.path.join(tables_out, f"{sample}_genes_cells.csv"))
    df_avg.to_csv(os.path.join(tables_out, f"{sample}_genes.csv"))

    # output dirs
    figs_out = os.path.join(figures_dir, out_folder)
    summary_out = os.path.join(figs_out, "summary_marked_genes")
    binary_out = os.path.join(figs_out, "binary_marked_genes")
    expr_out = os.path.join(figs_out, "binary_marked_genes")
    os.makedirs(summary_out, exist_ok=True)
    os.makedirs(binary_out, exist_ok=True)
    os.makedirs(expr_out, exist_ok=True)

    # --- UMAP embedding ---
    if "X_umap" not in adata.obsm:
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
        title=f"{sample}\nCell Type = boris_marker",
        palette=palette,
    )
    fig.set_size_inches(6, 6)
    fig.savefig(os.path.join(summary_out, f"celltype_{sample}_umap.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 3) Binary UMAPs: top 10 + always CTCFL
    if plot_yes:
        binary_genes = list(top_genes)
        ctcfl_matches = [g for g in all_genes_list if g.upper() == "CTCFL"]
        if ctcfl_matches:
            ctcfl_gene = ctcfl_matches[0]
            if ctcfl_gene not in binary_genes:
                binary_genes.append(ctcfl_gene)

        for g in binary_genes:
            if g not in df.columns:
                continue
            expr = df[g].values
            pos_mask = expr > 0
            n_pos, n_total = int(pos_mask.sum()), expr.shape[0]

            adata.obs[f"{g}_binary"] = pd.Categorical(
                np.where(pos_mask, "pos", "neg"),
                categories=["neg", "pos"]
            )

            pos_color = GENE_COLOR_MAP.get(g,"#d3d3d3")
            palette = {"neg": "#d3d3d3", "pos": pos_color if n_pos > 0 else "#d3d3d3"}

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

    return df_avg

# -------------------- run as script --------------------
if __name__ == "__main__":
    import scanpy as sc
    adata = sc.read_h5ad("/vf/users/scottaa/scrnaseq_pipeline/working_adata_tumors/T77.h5ad")
    export_gene_expression_and_umap(
        adata,
        sample="T77",
        tables_dir=TABLES_DIR_TUMOR,
        figures_dir=FIGURES_DIR_TUMOR,
        out_folder="boris_only",
        important_genes=ALL_GENES,  # or your gene list
        celltype_key="cell_type",
        plot_yes=True,
    )