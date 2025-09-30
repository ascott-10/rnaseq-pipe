#python_scripts/config.py

#### Libraries ###
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

#### Directories ###
BASE_DIR ="/vf/users/scottaa/scrnaseq_pipeline/"
REFERENCE_DIR = os.path.join(BASE_DIR,"reference_data")
MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR,"markers_df_fetal_gonad.csv")



#tumor files
RAW_DATA_DIR_TUMOR = os.path.join(BASE_DIR,"raw_data_tumors")
WORKING_DIR_TUMOR = os.path.join(BASE_DIR,"working_adata_tumors")
FIGURES_DIR_TUMOR = os.path.join(BASE_DIR,"results_tumors", "figures")
TABLES_DIR_TUMOR = os.path.join(BASE_DIR,"results_tumors", "tables")

#ovarian cancer files
RAW_DATA_DIR_OV = os.path.join(BASE_DIR,"raw_data_ovarian")
WORKING_DIR_OV = os.path.join(BASE_DIR,"working_adata_ovarian")
FIGURES_DIR_OV = os.path.join(BASE_DIR,"results_ovarian", "figures")
TABLES_DIR_OV = os.path.join(BASE_DIR,"results_ovarian", "tables")

#human fetal gonad files
RAW_DATA_DIR_FETAL = os.path.join(BASE_DIR,"raw_data_human_fetal_gonad")
WORKING_DIR_FETAL = os.path.join(BASE_DIR,"working_adata_human_fetal_gonad")
FIGURES_DIR_FETAL = os.path.join(BASE_DIR,"results_human_fetal_gonad", "figures")
TABLES_DIR_FETAL = os.path.join(BASE_DIR,"results_human_fetal_gonad", "tables")

#human fallopain tubes files
RAW_DATA_DIR_TUBES = os.path.join(BASE_DIR,"raw_data_tubes")
WORKING_DIR_TUBES = os.path.join(BASE_DIR,"working_adata_tubes")
FIGURES_DIR_TUBES = os.path.join(BASE_DIR,"results_human_tubes", "figures")
TABLES_DIR_TUBES = os.path.join(BASE_DIR,"results_human_tubes", "tables")

#human adult ovary files
RAW_DATA_DIR_ADULT_OVARY = os.path.join(BASE_DIR,"raw_data_adult_ovary")
WORKING_DIR_ADULT_OVARY = os.path.join(BASE_DIR,"working_adata_adult_ovary")
FIGURES_DIR_ADULT_OVARY = os.path.join(BASE_DIR,"results_adult_ovary", "figures")
TABLES_DIR_ADULT_OVARY = os.path.join(BASE_DIR,"results_adult_ovary", "tables")

#human embryo files
RAW_DATA_DIR_EMBRYOS = os.path.join(BASE_DIR,"raw_data_embryos")
WORKING_DIR_EMBRYOS = os.path.join(BASE_DIR,"working_adata_embryos")
FIGURES_DIR_EMBRYOS = os.path.join(BASE_DIR,"results_embryos", "figures")
TABLES_DIR_EMBRYOS = os.path.join(BASE_DIR,"results_embryos", "tables")

MENO_LIST = ["D1_3041", "D2_3061", "D3_3203", "D4_3295", "D5_3296", "D6_3369", "D7_3391"]
HEALTHY_LIST = ["donor_1", "donor_2", "donor_3", "donor_4"]

#### Parameters ####
RAW_FILE_TYPE_TUMOR = ".mtx" #Options: h5, mtx, csv
RAW_FILE_TYPE_OV = ".tsv" #Options: h5, mtx, csv, tsv
RAW_FILE_TYPE_FETAL = ".h5"
RAW_FILE_TYPE_TUBES = ".h5ad"
RAW_FILE_TYPE_ADOV = ".mtx"
RAW_FILE_TYPE_EMBRYOS = ".txt"
LOG2_CUTOFF = int(2)
LOG2_CUTOFF_2 = int(1)
PVAL_CUTOFF = float(0.001)
PVAL_CUTOFF_2 = float(0.05)

# Gene sets for analysis
IMPORTANT_GENES = ["CTCFL", "CTCF", "DPEP3", "BRCA2", "STAG3", "PRSS50", "MEIOC", "MMP10", "MAGEA12", "SAGE1", "PRAME", "WT1", "WISP2", "WNTB2", "ADRA2"]
GENES_CORRELATED = ['ACVR1C', 'ADD2', 'ARHGAP4', 'ARL2', 'C3ORF33', 'CCDC42', 'CCNJ', 'CLDN6', 'COL9A1', 'COQ10A', 'CREM', 'CT45A5', 'CYP26C1', 'DPEP3', 'F12', 'FOXA3', 'FXYD7', 'GCH1', 'GUSBP5', 'HIF3A', 'ING1', 'JADE3', 'KCNJ11', 'KLHL14', 'MAGEA4', 'MEIOC', 'NMUR1', 'NOX4', 'NTS', 'OBI1', 'PCCB', 'PCK1', 'PILRA', 'PLEKHG4', 'PRAME', 'PRSS21', 'RBM38', 'RSPO3', 'SCN8A', 'SELENOM', 'SFRP5', 'SLC38A2', 'SLC8A1', 'SLCO4C1', 'SPMIP5', 'TKTL1', 'TM7SF2', 'TMC7', 'TMEM44', 'TMTC4', 'TP53I3', 'TRIM58', 'UGT8', 'XXYLT1']
CTA_GENES = ["MAGEA4", "PRAME", "CT45A5", "MEIOC", "JADE3"]
METABOLIC_GENES = ["PCCB", "PCK1", "TKTL1", "COQ10A"]
SIGNALING_GENES = ["CREM", "FOXA3", "ING1", "HIF3A", "RSPO3", "SFRP5", "NMUR1", "NTS", "PILRA"]
TRANSPORT_GENES = ["ABCA3", "SLC38A2", "SLC8A1", "KCNJ11", "SCN8A"]
ECM_GENES = ["COL9A1", "CLDN6", "ADD2"]

ALL_GENES = list(set(
    IMPORTANT_GENES +
    GENES_CORRELATED +
    CTA_GENES +
    METABOLIC_GENES +
    SIGNALING_GENES +
    TRANSPORT_GENES +
    ECM_GENES
))

CTCFL_MARKERS = ['CREM', 'CLDN6', 'MAGEA4', 'DPEP3', 'TP53I3', 'WT1', 'CTCF', 'CTCFL', 'PRAME', 'ING1', 'RBM38', 'TKTL1', 'TM7SF2', 'RSPO3', 'XXYLT1', 'ARHGAP4', 'PRSS21', 'PCCB', 'HIF3A', 'F12', 'TMTC4', 'SLC8A1', 'NTS', 'COQ10A', 'CTCFL', 'COL9A1', 'WISP2', 'GCH1', 'PCK1', 'STAG3', 'BRCA2', 'FXYD7', 'MMP10', 'ADD2', 'KLHL14', 'FOXA3', 'CCNJ', 'JADE3', 'UGT8', 'MEIOC', 'ACVR1C', 'SAGE1']

ALL_GENES = list(set(
    IMPORTANT_GENES +
    GENES_CORRELATED +
    CTA_GENES +
    METABOLIC_GENES +
    SIGNALING_GENES +
    TRANSPORT_GENES +
    ECM_GENES
))


cells_df = pd.read_csv(MARKERS_FILE_PATH)
CELL_TYPES_LIST = list(cells_df.source.unique())



ALL_SAMPLES_TUMOR = ['T59', 'T76', 'T77', 'T89', 'T90']



ALL_SAMPLES_FETAL = ['mesenephros_M', 'mesenphros_F', 'mixed_G6_A', 'mixed_G6_B', 'ovary_1_G1', 'ovary_2_G2', 'ovary_G5_A', 'ovary_G5_B', 'testis_1_G3', 'testis_2_G4']

TUMOR_SAMPLES = ['T59', 'T76', 'T77', 'T89', 'T90']

OVARIAN_SAMPLES = ['OV-020-01-1A', 'OV-020-02-1A', 'OV-020-03-1A', 'OV-020-04-1A', 'OV-030-01-1A', 'OV-030-02-1A', 'OV-030-03-1A', 'OV-030-04-1A', 'OV-030-05-1A', 'OV-030-06-1A', 'OV-030-07-1A', 'OV-030-08-1A', 'OV-030-09-1A', 'OV-030-10-1A']
OVARIAN_DONOR_NAMES = ['BT1305', 'BT1306', 'BT1307', 'scrSOL001', 'scrSOL003', 'scrSOL004','scrSOL006', 'BT1303', 'BT1304','scrSOL007']
ALL_SAMPLES_TUBES = ['healthy_fallopian_raw_file', 'post_menopausal_raw_file']

ALL_SAMPLES_TUBE = ['healthy_fallopian_raw_file', 'post_menopausal_raw_file']

ALL_SAMPLES_ADULT_OVARY = ['middle-2', 'middle-3', 'middle-4', 'old-1', 'old-2', 'old-3', 'young-1', 'young-3', 'young-4']



EMBRYO_ALL_REPLICATES = ['F_10W_combined', 'F_10W_embryo1', 'F_12W_combined', 'F_12W_embryo1', 'F_14W_combined', 'F_14W_embryo1_1', 'F_14W_embryo1_2', 'F_14W_embryo1_3', 'F_18W_combined', 'F_18W_embryo1', 'F_18W_embryo2', 'F_20W_combined', 'F_20W_embryo1', 'F_20W_embryo2', 'F_23W_combined', 'F_23W_embryo1', 'F_23W_embryo2', 'F_24W_combined', 'F_24W_embryo1', 'F_24W_embryo2', 'F_26W_combined', 'F_26W_embryo1', 'F_5W_combined', 'F_5W_embryo1_and_2', 'F_7W_combined', 'F_7W_embryo1', 'F_8W_combined', 'F_8W_embryo1', 'M_10W_combined', 'M_10W_embryo1', 'M_10W_embryo2', 'M_12W_combined', 'M_12W_embryo1', 'M_19W_combined', 'M_19W_embryo1_101', 'M_19W_embryo1_24', 'M_19W_embryo1_26', 'M_19W_embryo2_102', 'M_19W_embryo2_103', 'M_19W_embryo2_104', 'M_20W_combined', 'M_20W_embryo1', 'M_21W_combined', 'M_21W_embryo1', 'M_21W_embryo2', 'M_21W_embryo3_10', 'M_21W_embryo3_10_2', 'M_21W_embryo3_17', 'M_25W_combined', 'M_25W_embryo1_101', 'M_25W_embryo1_102', 'M_25W_embryo1_103', 'M_25W_embryo1_104', 'M_25W_embryo1_105', 'M_25W_embryo1_106', 'M_25W_embryo1_107', 'M_25W_embryo1_24', 'M_25W_embryo1_26', 'M_4W_combined', 'M_4W_embryo1_and_F_11W_embryo1', 'M_9W_combined', 'M_9W_embryo1']

EMBRYO_SAMPLES = ['F_10W', 'F_12W', 'F_14W', 'F_18W', 'F_20W', 'F_23W', 'F_24W', 'F_26W', 'F_5W', 'F_7W', 'F_8W', 'M_10W', 'M_12W', 'M_19W', 'M_20W', 'M_21W', 'M_25W', 'M_4W', 'M_9W']
