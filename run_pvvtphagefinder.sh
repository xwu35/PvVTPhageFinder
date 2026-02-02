#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J run_PvVTPhageFinder
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xiaofen

eval "$(conda shell.bash hook)"
conda activate snakemake

python PvVTPhageFinder.py \
    --reads_dir /path/to/dir/containing/sequences \
    --sample_info /path/to/sample_info_table \
    --output_dir output_dir \
    --reference_genome /path/to/reference/genome/fasta/file 
