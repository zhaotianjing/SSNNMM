#! /bin/bash


pctAll=("0.99" "0.9" "0.7" "0.5")


for pct in "${pctAll[@]}"
do
	sbatch accuracy_nnmm_pblup_zw1.sbatch $pct
done