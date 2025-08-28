#!/bin/bash
#SBATCH --job-name=install_env
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/conda_%j.out
#SBATCH --error=logs/conda_%j.err

module load Anaconda3  

cd $HOME/rnaseq/
conda env create -f rnaseq_env.yaml
