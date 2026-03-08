# python_scripts_common/main.py

# Import libraries
import sys

# Custom
from config import get_project_config, PROJECT_CONFIG as project_config, BASE_DIR as default_base_dir
from main_process_data import run_process_data_pipeline

def main(project_name, base_dir = default_base_dir):

    

    if project_name not in project_config: 
        raise ValueError(f"Unknown project: {project_name}\n  Choose from {project_config.keys()}")
    
    # Otherwise set up directory paths
    paths = get_project_config(project_name, base_dir = default_base_dir)

    config = project_config[project_name]

    #Pas everything 
    run_process_data_pipeline(project_name = project_name, base_dir = base_dir, paths = paths, config = config)


if __name__ == "__main__":

    project_name = sys.argv[1] if len(sys.argv) > 1 else "mtab_tumors"
    base_dir = sys.argv[2] if len(sys.argv) > 2 else default_base_dir

    main(project_name, base_dir)