#!/bin/bash
#SBATCH --job-name=run_python       # Job name
#SBATCH --output=logs/output_%j.log         # Output file (%j will be replaced by job ID)
#SBATCH --error=errors/error_%j.log           # Error file
#SBATCH --time=01:00:00                # Time limit (hh:mm:ss)
#SBATCH --partition=norm            # Partition name
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=4G                       # Memory per node

# Load Python module (if required)
module load python/3.10

# Run your Python script
python my_script.py
