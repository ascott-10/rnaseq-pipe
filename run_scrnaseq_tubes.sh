#!/bin/bash
#SBATCH --job-name=tubes_scrnaseq
#SBATCH --output=logs/tubes_scrnaseq_%j.out
#SBATCH --error=errors/tubes_scrnaseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

cd /data/scottaa/scrnaseq_pub/projects/tubes

python python_scripts/main_process_raw_data.py

# To run (example):
# sbatch run_scrnaseq_tubes.sh
