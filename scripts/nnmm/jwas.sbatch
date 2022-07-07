#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --time=9-99:00:00
#SBATCH --partition=
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=



module load julia

srun julia --threads 10 /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/jwas.jl $1 $2 $3