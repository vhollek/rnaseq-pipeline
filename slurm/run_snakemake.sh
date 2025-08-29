#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err

module load Anaconda3  

snakemake \
    --jobs 100 \
    --use-conda \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster "sbatch --cpus-per-task={threads}" \
    --keep-going
