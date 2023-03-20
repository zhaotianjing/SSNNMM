#! /bin/bash

nRep=10  #10
methodAll=("ss") 
pctAll=("0.9" "0.7" "0.5" "0.3")   

cd /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000

for pct in "${pctAll[@]}"
do
	cd "$pct"

	for method in "${methodAll[@]}"
	do
		cd run_"$method"

		for rep in $( eval echo {1..$nRep} )  #rep1,...,rep20
		do
		    sbatch /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/ss/jwas_conventional_ss.sbatch $method $rep $pct
		done

		cd .. 
	done
	cd ..
done