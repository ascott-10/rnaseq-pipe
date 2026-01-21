# projects/E-MTAB-8559/python_scripts/config.py

import os

BASE_DIR = "/data/scottaa/scrnaseq_pub/"
REFERENCE_DIR = os.path.join(BASE_DIR, "reference_data")
PROJECT_NAME = "E-MTAB-8559"
PROJECT_DIR=os.path.join(BASE_DIR,"projects", PROJECT_NAME)

MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR, "markers_df.csv")
ALL_GENES_FILE_PATH = os.path.join(REFERENCE_DIR,  "ensembl_hgnc.csv")

RAW_DATA_DIR = os.path.join(PROJECT_DIR, "raw_data")
WORKING_ADATA_DIR = os.path.join(PROJECT_DIR, "working_adata")
RESULTS_DIR=os.path.join(PROJECT_DIR, "results")
FIGURES_DIR = os.path.join(PROJECT_DIR, "results", "figures")
TABLES_DIR = os.path.join(PROJECT_DIR, "results", "tables")

QC_DIR = os.path.join(FIGURES_DIR,"qc_violin_plot")
UMAP_DIR = os.path.join(FIGURES_DIR,"umap")

RAW_FILE_TYPE = ".mtx"

SAMPLE_METADATA = os.path.join(PROJECT_DIR, "sample_metadata.csv")
WHOLE_SAMPLE = [PROJECT_NAME]
SAMPLE_LIST = ["38b", "59", "74-1", "79"]


IMPORTANT_GENES = ['CTCFL',
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

for dir in [RAW_DATA_DIR, WORKING_ADATA_DIR,RESULTS_DIR,FIGURES_DIR,TABLES_DIR, QC_DIR, UMAP_DIR]:
    os.makedirs(dir, exist_ok=True)


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