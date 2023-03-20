#! /bin/bash

nRep=10
pctAll=("0.9" "0.7" "0.5" "0.3") 
method="centerzg_fixvar_noy"

cd /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000

for pct in "${pctAll[@]}"
do  
	cd "$pct"
	cd run_"$method"
	
	for rep in $( eval echo {1..$nRep} )  #rep1,...,rep20
	do
	    sbatch /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/"$method"/ssnnmm.sbatch $rep $pct $method
	done
	
	cd ..
	cd ..
done