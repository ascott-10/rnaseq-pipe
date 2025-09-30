#!/bin/bash
#SBATCH --output=logs/output_%j.log
#SBATCH --error=errors/error_%j.log
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G  

TYPE=$1
echo "Running sample type: $TYPE"

source myconda
conda activate rnaseq-pipe
python python_scripts/main.py $TYPE


# To run:

#sbatch --job-name=tumor run_python_sbatch.sh tumor_samples
#sbatch --job-name=ovarian run_python_sbatch.sh OVARIAN
#sbatch --job-name=fetal run_python_sbatch.sh fetal_samples

