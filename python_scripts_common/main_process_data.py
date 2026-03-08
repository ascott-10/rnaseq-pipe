# python_scripts_common/main_process_raw_data.py


#### Libraries ###

import os
import scanpy as sc

#### Custom ####

from config import *

from codes.process_raw_data import import_raw_data_mtx, run_qc, full_cell_annotation, make_umaps
#from codes.cell_annotation import make_universal_gene_list, init_color_maps, get_cell_type

### Process raw data ###

def run_process_data_pipeline(project_name, base_dir, paths, config):

    gene_map_file_path=ENSEMBL_TO_HGNC_MAP
    important_genes = IMPORTANT_GENES
    markers_file_path = MARKERS_FILE_PATH_GLOBAL
    cell_annotation_file_path = CELL_ANNOTATION_FILE_PATH

    cell_type_colors = CELL_TYPE_COLORS
    gene_color_map = GENE_COLOR_MAP

    raw_data_dir = paths["raw_data_dir"]
    working_adata_dir = paths["working_adata_dir"]
    umap_dir = paths["umap_dir"]
    qc_dir = paths["qc_dir"]

    raw_file_type = config["raw_counts_type"]
    project_type = config["project_type"]
    meta_df_filename = config["metadata_filename"]
    meta_df_file_path = os.path.join(paths["reference_data"], meta_df_filename)

    sample_list_dict = config["sample_list_dict"]
    if config["raw_data_dict"] is not None:
        raw_data_dict = config["raw_data_dict"]

    if project_type == "single_sample":
        # ["mtab-tumors"]
        adata_init_path = os.path.join(working_adata_dir, f"{project_name}_init.h5ad")
        adata_qc_path = os.path.join(working_adata_dir, f"{project_name}_after_qc.h5ad")
        adata_annotated_path = os.path.join(working_adata_dir, f"{project_name}_annotated.h5ad")

        if os.path.exists(adata_annotated_path):
            adata_annotated = sc.read_h5ad(adata_annotated_path)
                
        else:
            if os.path.exists(adata_qc_path):
                adata_after_qc = sc.read_h5ad(adata_qc_path)
            
            
            else: 



                if os.path.exists(adata_init_path):
                    adata_init = sc.read_h5ad(adata_init_path)
                else:
                    adata_init = import_raw_data_mtx(project_name, raw_data_dict,raw_data_dir,adata_init_path, meta_df_file_path,gene_map_file_path)

                adata_after_qc = run_qc(adata = adata_init, sample_id = project_name, adata_qc_path = adata_qc_path, qc_dir = qc_dir, important_genes = important_genes)


            adata_annotated = full_cell_annotation(sample = project_name, adata = adata_after_qc, adata_annotated_path = adata_annotated_path, important_genes = important_genes, markers_file_path = markers_file_path, cell_annotation_file_path = cell_annotation_file_path, cell_type_label_col="predicted_cell_type")
        
        
        
        
        make_umaps(sample = project_name, adata = adata_annotated, umap_dir = umap_dir, genes_list = important_genes, gene_color_map = gene_color_map, cell_type_colors = cell_type_colors)
           