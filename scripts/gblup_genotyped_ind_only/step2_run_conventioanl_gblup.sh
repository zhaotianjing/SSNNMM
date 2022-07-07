#! /bin/bash

nRep=10
methodAll=("gblup")
pctAll=("0.99" "0.9" "0.7" "0.5")


for pct in "${pctAll[@]}"
do
	cd "$pct"

	for method in "${methodAll[@]}"
	do
		cd run_"$method"

		for rep in $( eval echo {1..$nRep} )
		do
		    sbatch /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/jwas_conventional_gblup.sbatch $method $rep $pct
		done

		cd .. 
	done
	cd ..
done