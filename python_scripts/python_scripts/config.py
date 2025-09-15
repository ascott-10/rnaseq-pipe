#python_scripts/config.py

#### Libraries ###
import os

#### Directories ###
BASE_DIR ="/vf/users/scottaa/scrnaseq_pipeline/"
RAW_DATA_DIR = os.path.join(BASE_DIR,"raw_data")
WORKING_DIR = os.path.join(BASE_DIR,"working_adata")
REFERENCE_DIR = os.path.join(BASE_DIR,"reference_data")
FIGURES_DIR = os.path.join(BASE_DIR,"results", "figures")
TABLES_DIR = os.path.join(BASE_DIR,"results", "tables")

#### Parameters ####
RAW_FILE_TYPE = ".mtx" #Options: h5, mtx, csv
MARKERS_FILE_PATH = os.path.join(REFERENCE_DIR,"markers_df.csv")
ALL_SAMPLES = ['T77', 'T76', 'T90', 'T89', 'T59']
