#! /bin/bash


pctAll=("0.99" "0.9" "0.7" "0.5")


for pct in "${pctAll[@]}"
do
	sbatch /group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/calc_accuracy/accuracy_gblup.sbatch $pct
done