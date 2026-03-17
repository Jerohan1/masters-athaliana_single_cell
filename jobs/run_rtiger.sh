#!/bin/bash
#SBATCH --job-name=rtiger
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cores=4
#SBATCH --account athaliana_single_cell

source $(conda info --base)/etc/profile.d/conda.sh
conda activate rtiger

Rscript /faststorage/project/athaliana_single_cell/workspace/masters-athaliana_single_cell/scripts/rtiger.R
