#!/bin/bash

#this file is to make folders
nRep=20  #rep1,...,rep20
methodAll=("nnmm_pblup" "ss" "gblup")  
pctAll=("0.1" "0.3" "0.5" "0.7" "0.8" "0.9" "0.95" "0.99") 

for pct in "${pctAll[@]}"
do
	mkdir -p "$pct"
	cd "$pct"

	for method in "${methodAll[@]}"
	do
	    mkdir -p "$method"
		  cd "$method"
	    for j in $( eval echo {1..$nRep} ) 
	    do
	      mkdir -p rep"$j"
	    done
		cd ..

		mkdir -p run_"$method"
		mkdir -p accuracy_"$method"
	done
	mkdir -p accuracy_nnmm_pblup_zw1
	cd ..
	
done


mkdir -p pblup
cd pblup
mkdir -p run_pblup
mkdir -p accuracy_pblup
cd ..

mkdir -p gblup_full
cd gblup_full
mkdir -p run_gblup_full
mkdir -p accuracy_gblup_full
cd ..

mkdir -p nnmm_noomics
cd nnmm_noomics
mkdir -p run_nnmm_noomics
mkdir -p accuracy_nnmm_noomics
cd ..
