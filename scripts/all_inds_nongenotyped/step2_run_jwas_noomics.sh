#! /bin/bash

nRep=10


for rep in $( eval echo {1..$nRep} )
do
    sbatch jwas_noomics.sbatch $rep
done
