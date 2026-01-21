# projects/ovarian_cancer/python_scripts/main_process_raw_data.py


#### Libraries ###

import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc

import scipy.sparse as sp  
from scipy.sparse import csr_matrix
from scipy.io import mmread

import matplotlib.pyplot as plt

#### Custom ####

from config import *

from codes.process_raw_data import import_raw_data,run_qc,make_umaps
from codes.cell_annotation import make_universal_gene_list, init_color_maps, get_cell_type

### Make palettes ### 
CLEANED_GENES = make_universal_gene_list(dataset_file_path = ALL_GENES_FILE_PATH, cell_markers_file_path = MARKERS_FILE_PATH)
GENE_COLOR_MAP, CELL_TYPE_COLORS = init_color_maps(CLEANED_GENES,markers_file_path=MARKERS_FILE_PATH, adata=None)


### Process raw data ###


for project in WHOLE_SAMPLE:

    adata_init_path = os.path.join(WORKING_ADATA_DIR, f"{project}_init.h5ad")
    adata_qc_path = os.path.join(WORKING_ADATA_DIR, f"{project}_after_qc.h5ad")

    if os.path.exists(adata_qc_path):
        adata_qc = sc.read_h5ad(adata_qc_path)
    else:
        if os.path.exists(adata_init_path):
            adata_qc = run_qc(project, adata_init_path, adata_qc_path, important_genes=IMPORTANT_GENES)
        else:
            import_raw_data(project, adata_init_path, raw_file_type=RAW_FILE_TYPE, raw_data_dir=RAW_DATA_DIR, meta_file_path=SAMPLE_METADATA)
            adata_qc = run_qc(project, adata_init_path, adata_qc_path, important_genes=IMPORTANT_GENES)


    make_umaps(sample=project, adata_working_path = adata_qc_path, umap_dir = UMAP_DIR, important_genes=SIG_LIST)


for sample in SAMPLE_LIST:

    bdata_path = os.path.join(WORKING_ADATA_DIR, f"{sample}_after_qc.h5ad")

    if not os.path.exists(bdata_path):
        bdata = adata_qc[adata_qc.obs["sample"] == sample].copy()
        bdata.write(bdata_path)
    else:
        bdata = sc.read_h5ad(bdata_path)

    make_umaps(sample, adata_working_path=bdata_path, umap_dir=UMAP_DIR, important_genes=SIG_LIST)



