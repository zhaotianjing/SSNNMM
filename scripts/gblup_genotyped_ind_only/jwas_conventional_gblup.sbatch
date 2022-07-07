#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=10:00:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia

srun julia /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/jwas.jl $1 $2 $3
