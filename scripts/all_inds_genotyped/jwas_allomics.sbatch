#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --time=9-99:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia

srun julia --threads 10 jwas_allomics.jl $1