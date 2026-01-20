#!/bin/bash
#SBATCH --job-name=ov_car_scrnaseq
#SBATCH --output=logs/ov_car_scrnaseq_%j.out
#SBATCH --error=errors/ov_car_scrnaseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

source /data/scottaa/conda/etc/profile.d/conda.sh
conda activate rnaseq-pipe

cd /data/scottaa/scrnaseq_pub/projects/ovarian_cancer

python python_scripts/main_process_raw_data.py

# To run (example):
# sbatch run_scrnaseq_ovarian_cancer.sh
