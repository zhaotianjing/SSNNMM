#!/bin/bash -l
#SBATCH --job-name=mpitest
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2500
#SBATCH --time=90:00:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu


module load julia

srun julia /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/"$3"/ssnnmm.jl $1 $2 $3