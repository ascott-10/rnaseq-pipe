#main.py

#### Libraries ###

import sys


#### Custom ####

from config import *
import setup
import codes.process_raw_data


#### Custom ####
import setup
import codes.process_raw_data
import codes.gene_expression
from config import MARKERS_FILE_PATH
from setup import GENE_COLOR_MAP, CELL_TYPE_COLORS   # if you split colors into colors.py
from codes.process_raw_data import make_adata              # wherever prep_adata/make_adata live
from setup import get_cell_type, get_raw_gene_counts # your helper functions

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main.py [setup | process_data <project_name>]")
        sys.exit(1)

    command = sys.argv[1]

    if command == "setup":
        
        setup.directory_setup()

    elif command == "process_data":
        if len(sys.argv) < 3:
            print("Usage: python main.py process_data <project_name>")
            sys.exit(1)
            
        GENE_COLOR_MAP, CELL_TYPE_COLORS = setup.init_color_maps(markers_file_path=MARKERS_FILE_PATH)
        current_project = sys.argv[2]
        raw_data_dir, working_adata_dir, figures_dir, tables_dir, raw_file_type, markers_file_path = setup.set_active_project(current_project)

        ALL_SAMPLES = codes.process_raw_data.detect_multiple_samples(raw_data_dir, current_project)
        for sample in ALL_SAMPLES:
            print(f"Processing sample {sample}")

            # load raw
            adata_init, adata_init_path = codes.process_raw_data.import_raw_data(
                sample, raw_data_dir, raw_file_type=raw_file_type
            )
            adata_working_path = os.path.join(working_adata_dir, f"{sample}.h5ad")

            # process
            adata = make_adata(
                adata_init_path,
                adata_working_path,
                sample,
                figures_dir,
                tables_dir,
                forced_genes=["CTCF", "CTCFL"],
                gene_color_map=GENE_COLOR_MAP,
                cell_type_colors=CELL_TYPE_COLORS,
                markers_file_path=MARKERS_FILE_PATH,
                get_cell_type_fn=get_cell_type,
                get_raw_gene_counts_fn=get_raw_gene_counts,
            )

            df_important_genes = codes.gene_expression.export_gene_expression_and_umap(adata,
                                                                                       sample,
                                                                                       tables_dir,
                                                                                       figures_dir,
                                                                                       out_folder = "test_1",
                                                                                       important_genes=["PLAC1"],
                                                                                       celltype_key="cell_type", 
                                                                                       plot_yes=True)
    else:
        print(f"Unknown command: {command}")




### Example Usage ###
#python main.py setup 



# To run:

#sbatch --job-name=tumor run_python_sbatch.sh tumor_samples
#sbatch --job-name=ovarian run_python_sbatch.sh OVARIAN
#sbatch --job-name=fetal run_python_sbatch.sh fetal_samples
#sbatch --job-name=tubes run_python_sbatch.sh tube_samples
#sbatch --job-name=adult run_python_sbatch.sh adult_ovary_samples
#sbatch --job-name=embryo run_python_sbatch.sh embryo_samples
