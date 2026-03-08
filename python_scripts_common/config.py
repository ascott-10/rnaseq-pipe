# python_scripts_common/config.py

import os

BASE_DIR = "/home/ascott10/documents/projects/rnaseq-pipe"

REFERENCE_DIR = os.path.join(BASE_DIR, "reference_data_common")
MARKERS_FILE_PATH_GLOBAL = os.path.join(REFERENCE_DIR, "CellMarker_markers_df.csv")
CELL_ANNOTATION_FILE_PATH = os.path.join(REFERENCE_DIR, "CellMarker_cell_type_annotation.csv")
ENSEMBL_TO_HGNC_MAP= os.path.join(REFERENCE_DIR,  "ensembl_to_hgnc.csv")
IMPORTANT_GENES = ["CTCF", "CTCFL", "DPEP3"]

GENE_DEFAULT_PALETTE = {
    "CTCFL": "#04531c",
    "CTCF": "#4b11d3",
    "DPEP3": "#d31162",
}

GENE_COLOR_MAP = GENE_DEFAULT_PALETTE.copy()
CELL_TYPE_COLORS = {"None": "#d3d3d3"}

PROJECT_CONFIG = {

    "mtab_tumors":
    
    {
        "markers_file_path": MARKERS_FILE_PATH_GLOBAL,
        "raw_counts_type": ".mtx",
        "metadata_filename" : "sample_metadata.csv",
        "project_type": "single_sample",            #single_sample (already concat). multi_sample (not pre-concat)
        "sample_list": ["38b", "59", "74-1", "79"], 
        "sample_list_dict": {"mtab_tumors": ["38b", "59", "74-1", "79"]},
         "raw_data_dict": { 
            "mtx_filename": "E-MTAB-8559.aggregated_filtered_counts.mtx",
            "cells_filename": "E-MTAB-8559.aggregated_filtered_counts.mtx_cols",
            "genes_filename" :"E-MTAB-8559.aggregated_filtered_counts.mtx_rows",
        },
    }

}


def get_project_config(project_name, base_dir = BASE_DIR):

    dirs = {
        "raw_data_dir": "raw_data",
        "working_adata_dir": "working_adata",
        "results_dir": "results",
        "figures_dir": "results/figures",
        "tables_dir": "results/tables",
        "qc_dir": "results/figures/qc",
        "umap_dir": "results/figures/umap",
        "reference_data": "reference_data"
    }

    paths = {}

    for dir_name, relative_path in dirs.items():
        full_path = os.path.join(base_dir, project_name, relative_path)
        os.makedirs(full_path, exist_ok=True)
        paths[dir_name] = full_path

    return paths


IMPORTANT_GENES_2 = ['CTCFL',
 'DPEP3',
 'SH3KBP1',
 'PCK1',
 'CKAP2',
 'PRSS50',
 'CEP55',
 'LAS1L',
 'LRP2',
 'PLEKHG4',
 'INTS6',
 'PHF16',
 'DNM1',
 'SEMA5B',
 'PARG',
 'MMP10',
 'TRIM17',
 'RPGR',
 'TMEM44',
 'PDE3A',
 'PDE11A',
 'ADRA2B',
 'GCH1',
 'RTKN2',
 'DBF4B',
 'STAG3',
 'WDR52',
 'ACTL8',
 'CCDC114',
 'AGAP6',
 'FXYD7',
 'SSX4',
 'BRCA2',
 'MYBL1',
 'CEP152',
 'EP400NL']




SIG_LIST =['ACTL8',
 'ATP2A1',
 'AXIN2',
 'BMP2',
 'BMP4',
 'BMP7',
 'BREA2',
 'CGB2',
 'CLDN3',
 'CLDN4',
 'CTNNB1',
 'DKK1',
 'DPEP3',
 'EHF',
 'ELF3',
 'EMX2',
 'EPCAM',
 'FOXA1',
 'FOXA2',
 'FSD1',
 'GATA3',
 'GPR3',
 'HOXA10',
 'HOXA11',
 'HOXA12',
 'HOXA13',
 'HOXA9',
 'INHBE',
 'KCTD19',
 'KRT18',
 'KRT19',
 'KRT8',
 'LCNL1',
 'LEF1',
 'LHX1',
 'LOC541473',
 'LRP2',
 'MMP10',
 'MUC1',
 'NCRNA00202',
 'PAX8',
 'PAX8 ',
 'PCK1',
 'PDE11A',
 'PHF16',
 'POPDC3',
 'PRSS50',
 'RAG1',
 'RTKN2',
 'SEMA5B',
 'SMAD2',
 'SMAD3',
 'SMAD4',
 'SSX4',
 'STAG3',
 'TCAM1P',
 'TCF7L2',
 'TCP10L',
 'TGFB1',
 'TGFB2',
 'TMC7',
 'TRIM17',
 'WNT4',
 'WNT7',
 'WNT7A',
 'WT1']