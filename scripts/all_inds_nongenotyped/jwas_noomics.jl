# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

rep          = parse(Int, ARGS[1])  #rep1,...,rep20
chain_length = 5000
ntest        = 0 #all inds are used for training, prediction accuracy are based on all individuals

@show rep,chain_length,ntest

mainpath="/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/no_omics/"
cd(mainpath)

################################################
# NN-PBLUP
################################################
############ READ DATA ##########
#read phenotype
phenofile="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/all_no_missing/y.rep$rep.csv"
phenotypes = CSV.read(phenofile,DataFrame);

#read pedigree
pedfile    = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/pedigree.txt"
pedigree   = get_pedigree(pedfile,separator=",",header=true);

# build model
model_equation = "y = intercept + ID";
model = build_model(model_equation;
                    num_hidden_nodes=5,
                    nonlinear_function="linear");
set_random(model,"ID",pedigree); 

out   = runMCMC(model,phenotypes,chain_length=chain_length,printout_model_info=false,estimate_variance=true,output_folder="nnmm.rep$rep");



################################################
# conventional PBLUP
################################################
############ READ DATA ##########
#read pedigree
pedfile    = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/pedigree.txt"
pedigree   = get_pedigree(pedfile,separator=",",header=true);

#read phenotype
phenofile="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/all_no_missing/y.rep$rep.csv"
phenotypes = CSV.read(phenofile,DataFrame);
phenotypes[!,:ID] = string.(phenotypes[!,:ID]);

# build model
model_equation = "y = intercept + ID";
model = build_model(model_equation);
set_random(model,"ID",pedigree); 

out   = runMCMC(model,phenotypes,chain_length=chain_length,output_folder="pblup.rep$rep");



