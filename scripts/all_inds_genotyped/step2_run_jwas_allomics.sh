#! /bin/bash

nRep=10  #3


for rep in $( eval echo {1..$nRep} )
do
    sbatch jwas_allomics.sbatch $rep
done
