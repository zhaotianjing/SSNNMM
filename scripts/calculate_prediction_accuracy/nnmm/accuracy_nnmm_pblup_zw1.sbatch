#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=00:20:00
#SBATCH --partition=
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=



module load julia

srun julia accuracy_nnmm_pblup_zw1.jl $1