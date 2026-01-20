#!/bin/bash
#SBATCH --job-name=fetal_scrnaseq
#SBATCH --output=logs/fetal_scrnaseq_%j.out
#SBATCH --error=errors/fetal_scrnaseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

cd /data/scottaa/scrnaseq_pub/projects/human_fetal_gonad

python python_scripts/main_process_raw_data.py

# To run (example):
# sbatch run_scrnaseq_fetal_gonad.sh
