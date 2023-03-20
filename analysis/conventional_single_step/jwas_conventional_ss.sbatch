#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=9G
#SBATCH --time=10:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu



module load julia/1.7.0

srun julia /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/ss/jwas.jl $1 $2 $3
