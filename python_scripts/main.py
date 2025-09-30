#### Libraries ###
import os
import gzip
import shutil

import pandas as pd
import scanpy as sc
import numpy as np
import decoupler as dc

import matplotlib.pyplot as plt
import matplotlib as mpl

import scipy.sparse as sp

# Scanpy / matplotlib defaults
 # bigger canvas
sc.settings.autoshow = False                                  # don't pop/show


# Local imports
from config import *


from codes.process_raw_data import (
    detect_multiple_samples, detect_replicates, combine_replicates,
    detect_compressed_files,
    import_raw_data,
    make_adata, make_bdata
)
from codes.utils import subset_adata, get_cell_type, add_list_to_config, show_examples, check_sources, get_raw_gene_counts, init_color_maps, GENE_COLOR_MAP



GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CTCFL_MARKERS + ALL_GENES)



from codes.sc_enrichment import (
    tf_scoring,
    pathway_scoring,
    hallmark_genesets,
    ctcfl_pathway_matrices,
    ctcfl_geneset_matrices,
)
from codes.psuedobulk import (
    setup_pseudobulk, singlecell_de_volcano,
    map_metadata,
    pseudobulk_dea,
    make_deseq_ready,
)
from codes.gene_expression import (
    export_all_expression,
    export_gene_split,
    plot_gene_expression,export_gene_expression_and_umap, find_coexpressed_genes
)


## Functions ##
def embryo_samples(raw_data_dir, raw_file_type, all_replicates_list_name, all_samples_list_name, working_dir, figures_dir, tables_dir, log2_cutoff, pval_cutoff):
    detect_compressed_files(raw_data_dir)
    EMBRYO_ALL_REPLICATES = detect_multiple_samples(raw_data_dir, list_name = "EMBRYO_ALL_REPLICATES")
    for replicate in EMBRYO_ALL_REPLICATES:
        
        if replicate.endswith("_combined"):  # skip combined pseudo-replicates
            continue

    

        print(f"Processing embryo replicate {replicate}")
        adata_init, adata_init_path = import_raw_data(replicate,raw_file_type,raw_data_dir) # Load raw data for this sample
                
    
    
    EMBRYO_SAMPLES = detect_replicates(raw_data_dir, base_list_name = "EMBRYO_SAMPLES")

    for sample in EMBRYO_SAMPLES:
        if sample.endswith("_combined"):  # skip combined pseudo-replicates
            continue
        adata_working_path = os.path.join(working_dir, f"{sample}.h5ad")

        if os.path.exists(adata_working_path):
            print(f"[{sample}] combined AnnData already exists")
            

        reps = [r for r in EMBRYO_ALL_REPLICATES if r.startswith(sample + "_")]

        if reps:

            print(f"[{sample}] combining {len(reps)} replicates…")
            adata_init = combine_replicates(sample, raw_data_dir)

            # save merged raw AnnData so make_adata works unchanged
            adata_init_path = os.path.join(raw_data_dir, f"{sample}_combined.h5ad")
            adata_init.write(adata_init_path)
        else:
            print(f"[{sample}] single file mode")
            adata_init_path = os.path.join(raw_data_dir, f"{sample}.h5ad")

        adata = make_adata(
                            adata_init_path,
                            adata_working_path,
                            sample,
                            figures_dir,
                            tables_dir,
                            # --- QC params with relaxed parameters ---
                            min_genes=1,
                            min_cells=1,
                            max_genes=6000,
                            max_mt=25,
                            # Outlier MAD thresholds
                            
                            
                            mt_prefix="MT-",
                            ribo_prefix=("RPS", "RPL"),
)

        #adata = get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)
        adata.write(adata_working_path)
        
        print(sample, adata.shape)
        print(adata)
        print(adata.obs.head())
        df_important_genes = export_gene_expression_and_umap(adata,sample,tables_dir,figures_dir,out_folder = "test_7",important_genes=ALL_GENES,celltype_key="cell_type", plot_yes=True)
        
        find_coexpressed_genes(adata,sample,target_gene = "CTCFL",tables_dir = tables_dir,figures_dir = figures_dir,out_folder = "test_7")


def ov_samples(raw_data_dir, raw_file_type, list_name, working_dir, figures_dir, tables_dir, log2_cutoff, pval_cutoff):
    detect_compressed_files(raw_data_dir)
    OVARIAN_SAMPLES = detect_multiple_samples(raw_data_dir, list_name)
    for sample in OVARIAN_SAMPLES:
    

        print(f"Processing ovarian sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample,raw_file_type,raw_data_dir) # Load raw data for this sample
        adata_working_path = os.path.join(working_dir, f"{sample}.h5ad")                   # Path for working copy
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir, tables_dir)    # Make working AnnData
        
        adata = get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)
        adata.write(adata_working_path)
        adata.obs["sample"] = sample  # ov-only customization
        
        df_important_genes = export_gene_expression_and_umap(adata,sample,tables_dir,figures_dir,out_folder = "test_7",important_genes=ALL_GENES,celltype_key="cell_type", plot_yes=True)
        
        find_coexpressed_genes(adata,sample,target_gene = "CTCFL",tables_dir = tables_dir,figures_dir = figures_dir,out_folder = "test_7")


        # pseudobulk + DEA
        #adata = setup_pseudobulk(adata, sample, tables_dir=tables_dir, figures_dir=figures_dir)
        #adata_counts = make_deseq_ready(adata, sample,min_gene_total=5, min_sample_total=2000, n_top_genes=None)
        #pseudobulk_dea(adata_counts, sample,log2_cutoff=log2_cutoff, pval_cutoff=pval_cutoff,tables_dir=tables_dir, figures_dir=figures_dir)

        #singlecell_de_volcano(adata, sample, figures_dir, tables_dir, groupby="CTCFL_expr",method="wilcoxon", log2_cutoff=1, pval_cutoff=0.05)



def tumor_samples(raw_data_dir, raw_file_type, list_name, working_dir, figures_dir, tables_dir, log2_cutoff, pval_cutoff):
    """Process tUMOR samples (patient tumors)."""
    TUMOR_SAMPLES = detect_multiple_samples(raw_data_dir, list_name)
    for sample in TUMOR_SAMPLES:
        print(f"Processing ovarian sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample,raw_file_type,raw_data_dir) # Load raw data for this sample
        adata_working_path = os.path.join(working_dir, f"{sample}.h5ad")                   # Path for working copy
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir, tables_dir)    # Make working AnnData
        
        adata = get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)
        adata.write(adata_working_path)
        # tumor-only customization: add clinical metadata
        adata = map_metadata(
            adata, sample, adata_working_path,
            meta_columns=("tumor_stage", "chemotherapy_response"),
            reference_dir=REFERENCE_DIR)

        df_important_genes = export_gene_expression_and_umap(adata,sample,tables_dir,figures_dir,out_folder = "test_7",important_genes=ALL_GENES,celltype_key="cell_type", plot_yes=True)
        
        find_coexpressed_genes(adata,sample,target_gene = "CTCFL",tables_dir = tables_dir,figures_dir = figures_dir,out_folder = "test_7")

        # pseudobulk + DEA
        #adata = setup_pseudobulk(adata, sample, tables_dir=tables_dir, figures_dir=figures_dir)
        #adata_counts = make_deseq_ready(adata, sample,min_gene_total=5, min_sample_total=2000, n_top_genes=None)
        #pseudobulk_dea(adata_counts, sample,log2_cutoff=log2_cutoff, pval_cutoff=pval_cutoff,tables_dir=tables_dir, figures_dir=figures_dir)

        #singlecell_de_volcano(adata, sample, figures_dir, tables_dir, groupby="CTCFL_expr",method="wilcoxon", log2_cutoff=1, pval_cutoff=0.05)
     

def fetal_samples(raw_data_dir, raw_file_type, list_name, working_dir, figures_dir, tables_dir, log2_cutoff, pval_cutoff):
    ALL_SAMPLES_FETAL = detect_multiple_samples(raw_data_dir, list_name)
    for sample in ALL_SAMPLES_FETAL:
   

   
    
        print(f"Processing FETAL sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample,raw_file_type,raw_data_dir) # Load raw data for this sample
        adata_working_path = os.path.join(working_dir, f"{sample}.h5ad")                   # Path for working copy
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir, tables_dir)    # Make working AnnData
        adata.obs["sample"] = sample  
        adata = get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)
        adata.write(adata_working_path)
        df_important_genes = export_gene_expression_and_umap(adata,sample,tables_dir,figures_dir,out_folder = "test_7",important_genes=ALL_GENES,celltype_key="cell_type", plot_yes=True)
        
        find_coexpressed_genes(adata,sample,target_gene = "CTCFL",tables_dir = tables_dir,figures_dir = figures_dir,out_folder = "test_7")


        # pseudobulk + DEA
        #adata = setup_pseudobulk(adata, sample, tables_dir=tables_dir, figures_dir=figures_dir)
        #adata_counts = make_deseq_ready(adata, sample,min_gene_total=5, min_sample_total=2000, n_top_genes=None)
        #pseudobulk_dea(adata_counts, sample,log2_cutoff=log2_cutoff, pval_cutoff=pval_cutoff,tables_dir=tables_dir, figures_dir=figures_dir)

        #singlecell_de_volcano(adata, sample, figures_dir, tables_dir, groupby="CTCFL_expr",method="wilcoxon", log2_cutoff=1, pval_cutoff=0.05)

def tube_samples(raw_data_dir, raw_file_type, list_name, working_dir, figures_dir, tables_dir, log2_cutoff, pval_cutoff):
    
    """Process Fallopian tube samples (Fallopian tubes)."""
    TUBE_SAMPLES = detect_multiple_samples(raw_data_dir, list_name)
    for sample in TUBE_SAMPLES:
        print(f"Processing Fallopian tube sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample,raw_file_type,raw_data_dir) # Load raw data for this sample
        adata_working_path = os.path.join(working_dir, f"{sample}.h5ad")                   # Path for working copy
        
       

        # donors actually present in this sample

        # Case-insensitive search for CTCFL in var_names
        matches = [g for g in adata_init.var_names if "CTCFL" in g.upper()]
        print("Matches for 'CTCFL':", matches)
        print("Global CTCFL max in adata_init:", adata_init[:, "CTCFL"].X.max())
        donors_in_sample = set(adata_init.obs["donor_id"].unique())

        for donor in HEALTHY_LIST + MENO_LIST:
            if donor not in donors_in_sample:
                print(f"⚠️ Skipping donor {donor}: not found in sample {sample}")
                continue
                
            print(f"Processing Fallopian tube donor {donor}")
            bdata_working_path = os.path.join(working_dir, f"{donor}.h5ad")
            bdata = make_bdata(donor, adata_init_path, adata_working_path, bdata_working_path, figures_dir, tables_dir)
            

            # Case-insensitive search for CTCFL in var_names
            # Same search inside one donor bdata
            matches_bdata = [g for g in bdata.var_names if "CTCFL" in g.upper() or "ENSG00000163631" in g.upper()]
            print("Matches in bdata:", donor, matches_bdata)

            print("Max in counts layer:", bdata.layers["counts"].max())

            bdata.write(bdata_working_path)
            df_important_genes = export_gene_expression_and_umap(bdata,donor,tables_dir,figures_dir,out_folder = "test_7",important_genes=ALL_GENES,celltype_key="cell_type", plot_yes=True)
        
            find_coexpressed_genes(bdata,donor,target_gene = "CTCFL",tables_dir = tables_dir,figures_dir = figures_dir,out_folder = "test_7")

        #singlecell_de_volcano(adata, sample, figures_dir, tables_dir, groupby="CTCFL_expr",method="wilcoxon", log2_cutoff=1, pval_cutoff=0.05)
 
def adult_ovary_samples(raw_data_dir, raw_file_type, list_name, working_dir, figures_dir, tables_dir, log2_cutoff, pval_cutoff):
    detect_compressed_files(raw_data_dir)
    ADULT_OVARY_SAMPLES = detect_multiple_samples(raw_data_dir, list_name)
    for sample in ADULT_OVARY_SAMPLES:
        print(f"Processing adult ovary sample {sample}")
        adata_init, adata_init_path = import_raw_data(sample,raw_file_type,raw_data_dir) # Load raw data for this sample
        adata_working_path = os.path.join(working_dir, f"{sample}.h5ad")                   # Path for working copy
        adata = make_adata(adata_init_path, adata_working_path, sample, figures_dir, tables_dir)    # Make working AnnData
        
        adata = get_cell_type(adata, markers_file_path = MARKERS_FILE_PATH)
        adata.write(adata_working_path)
        adata.obs["sample"] = sample  # ov-only customization
        
        df_important_genes = export_gene_expression_and_umap(adata,sample,tables_dir,figures_dir,out_folder = "test_7",important_genes=["PLAC1"],celltype_key="cell_type", plot_yes=True)
        find_coexpressed_genes(adata,sample,target_gene = "CTCFL",tables_dir = tables_dir,figures_dir = figures_dir,out_folder = "test_7")


    

## Entry Point ##

import sys

def main(sample_type):
    if sample_type == "tumor_samples":
        tumor_samples(
            RAW_DATA_DIR_TUMOR, RAW_FILE_TYPE_TUMOR, "TUMOR_SAMPLES",
            WORKING_DIR_TUMOR, FIGURES_DIR_TUMOR, TABLES_DIR_TUMOR,
            LOG2_CUTOFF_2, PVAL_CUTOFF_2
        )
    elif sample_type == "OVARIAN":
        ov_samples(
            RAW_DATA_DIR_OV, RAW_FILE_TYPE_OV, "OVARIAN_SAMPLES",
            WORKING_DIR_OV, FIGURES_DIR_OV, TABLES_DIR_OV,
            LOG2_CUTOFF_2, PVAL_CUTOFF_2
        )

    elif sample_type == "embryo_samples":
        embryo_samples(
            RAW_DATA_DIR_EMBRYOS, RAW_FILE_TYPE_EMBRYOS, "EMBRYO_ALL_REPLICATES","EMBRYO_SAMPLES",
            WORKING_DIR_EMBRYOS, FIGURES_DIR_EMBRYOS, TABLES_DIR_EMBRYOS,
            LOG2_CUTOFF_2, PVAL_CUTOFF_2
        )
    elif sample_type == "fetal_samples":
        fetal_samples(
            RAW_DATA_DIR_FETAL, RAW_FILE_TYPE_FETAL, "ALL_SAMPLES_FETAL",
            WORKING_DIR_FETAL, FIGURES_DIR_FETAL, TABLES_DIR_FETAL,
            LOG2_CUTOFF_2, PVAL_CUTOFF_2
        )
    elif sample_type == "tube_samples":
        tube_samples(RAW_DATA_DIR_TUBES, RAW_FILE_TYPE_TUBES, "ALL_SAMPLES_TUBE",
            WORKING_DIR_TUBES, FIGURES_DIR_TUBES, TABLES_DIR_TUBES,
            LOG2_CUTOFF_2, PVAL_CUTOFF_2
        )

    elif sample_type == "adult_ovary_samples":
        
        
        adult_ovary_samples(RAW_DATA_DIR_ADULT_OVARY, RAW_FILE_TYPE_ADOV, "ALL_SAMPLES_ADULT_OVARY",
            WORKING_DIR_ADULT_OVARY, FIGURES_DIR_ADULT_OVARY, TABLES_DIR_ADULT_OVARY,
            LOG2_CUTOFF_2, PVAL_CUTOFF_2        )
    else:
        raise ValueError("sample_type must be one of: tumor_samples, OVARIAN, fetal_samples, tube_samples", "adult_ovary_samples")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError(
            "Usage: python main.py <sample_type>\n"
            
        )
    sample_type = sys.argv[1]
    main(sample_type)

# To run:

#sbatch --job-name=tumor run_python_sbatch.sh tumor_samples
#sbatch --job-name=ovarian run_python_sbatch.sh OVARIAN
#sbatch --job-name=fetal run_python_sbatch.sh fetal_samples
#sbatch --job-name=tubes run_python_sbatch.sh tube_samples
#sbatch --job-name=adult run_python_sbatch.sh adult_ovary_samples
#sbatch --job-name=embryo run_python_sbatch.sh embryo_samples
