#python_scripts/config.py

#### Libraries ###
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

#### Directories ###

import os

BASE_DIR = "/vf/users/scottaa/scrnaseq_pub/"
REFERENCE_DIR = os.path.join(BASE_DIR, "reference_data")
MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR, "markers_df.csv")

# User inputs
PROJECTS = {
    "tumor_samples": {
        "raw_file_type": ".mtx",
        "markers_file_path": os.path.join(REFERENCE_DIR, "markers_df.csv"),
    },
    "fetal_gonad": {
        "raw_file_type": ".h5",
        "markers_file_path": os.path.join(REFERENCE_DIR, "markers_df.csv")
    },
}



# template names
RAW_DATA_DIR     = "raw_data"
WORKING_DATA_DIR = "working_adata"
FIGURES_DIR      = os.path.join("results", "figures")
TABLES_DIR       = os.path.join("results", "tables")




IMPORTANT_GENES = ["CTCFL", "CTCF", "DPEP3", "BRCA2", "STAG3", "PRSS50", "MEIOC", "MMP10", "MAGEA12", "SAGE1", "PRAME", "WT1", "WISP2", "WNTB2", "ADRA2"]
GENES_CORRELATED = ['ACVR1C', 'ADD2', 'ARHGAP4', 'ARL2', 'C3ORF33', 'CCDC42', 'CCNJ', 'CLDN6', 'COL9A1', 'COQ10A', 'CREM', 'CT45A5', 'CYP26C1', 'DPEP3', 'F12', 'FOXA3', 'FXYD7', 'GCH1', 'GUSBP5', 'HIF3A', 'ING1', 'JADE3', 'KCNJ11', 'KLHL14', 'MAGEA4', 'MEIOC', 'NMUR1', 'NOX4', 'NTS', 'OBI1', 'PCCB', 'PCK1', 'PILRA', 'PLEKHG4', 'PRAME', 'PRSS21', 'RBM38', 'RSPO3', 'SCN8A', 'SELENOM', 'SFRP5', 'SLC38A2', 'SLC8A1', 'SLCO4C1', 'SPMIP5', 'TKTL1', 'TM7SF2', 'TMC7', 'TMEM44', 'TMTC4', 'TP53I3', 'TRIM58', 'UGT8', 'XXYLT1']
CTA_GENES = ["MAGEA4", "PRAME", "CT45A5", "MEIOC", "JADE3"]
METABOLIC_GENES = ["PCCB", "PCK1", "TKTL1", "COQ10A"]
SIGNALING_GENES = ["CREM", "FOXA3", "ING1", "HIF3A", "RSPO3", "SFRP5", "NMUR1", "NTS", "PILRA"]
TRANSPORT_GENES = ["ABCA3", "SLC38A2", "SLC8A1", "KCNJ11", "SCN8A"]
ECM_GENES = ["COL9A1", "CLDN6", "ADD2"]

#https://www.oatext.com/pdf/ICST-3-192.pdf
#https://www.nature.com/articles/s42003-024-06689-2
OTHER_GENES = ["PLAC1", "BAP1"]

CTCFL_MARKERS = ['CREM', 'CLDN6', 'MAGEA4', 'DPEP3', 'TP53I3', 'WT1', 'CTCF', 'CTCFL', 'PRAME', 'ING1', 'RBM38', 'TKTL1', 'TM7SF2', 'RSPO3', 'XXYLT1', 'ARHGAP4', 'PRSS21', 'PCCB', 'HIF3A', 'F12', 'TMTC4', 'SLC8A1', 'NTS', 'COQ10A', 'CTCFL', 'COL9A1', 'WISP2', 'GCH1', 'PCK1', 'STAG3', 'BRCA2', 'FXYD7', 'MMP10', 'ADD2', 'KLHL14', 'FOXA3', 'CCNJ', 'JADE3', 'UGT8', 'MEIOC', 'ACVR1C', 'SAGE1']

ALL_GENES = list(set(
    IMPORTANT_GENES +
    GENES_CORRELATED +
    CTCFL_MARKERS
))
tumor_samples = ['T59', 'T76', 'T77', 'T89', 'T90']

fetal_gonad = ['mesenephros_M', 'mesenphros_F', 'mixed_G6_A', 'mixed_G6_B', 'ovary_1_G1', 'ovary_2_G2', 'ovary_G5_A', 'ovary_G5_B', 'testis_1_G3', 'testis_2_G4']
