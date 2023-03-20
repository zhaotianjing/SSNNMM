# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

using Pkg
@show Pkg.status() #check JWAS version

println("-----02172023-----")
method            = ARGS[1]  # ("nnmm_pblup" "ss" "nnmm_gblup")
rep               = parse(Int, ARGS[2])  #rep1,...,rep20
pct               = parse(Float64, ARGS[3])
chain_length      = 5000
ntest=100 #number of testing individuals

@show method,rep,pct,chain_length,ntest

mainpath="/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/"
cd(mainpath*"$pct/$method/rep$rep/")


#conventional single-step
if method=="ss"
    println("-------conventional single-step analysis-------")
    ############ READ DATA ##########
    #read partial genotype
    genofile   ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p10000/single_step_partial_genotype/partialgenotype.pct$pct.rep$rep.csv"
    genotypes  = get_genotypes(genofile,separator=',',method="RR-BLUP",quality_control=false);

    #read pedigree
    pedfile    = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p10000/pedigree.txt"
    pedigree   = get_pedigree(pedfile,separator=",",header=true);

    #read phenotype
    phenofile  ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p10000/y.csv"
    phenotypes = CSV.read(phenofile,DataFrame);
    phenotypes[!,:ID] = string.(phenotypes[!,:ID]);

    phenotypes[!,:y]= convert(Vector{Union{Missing,Float64}}, phenotypes[!,:y]); #convert the data type to allow missing data
    phenotypes[end-ntest+1:end, :y] .= missing;  #missing phenotype for individuals in testing dataset (last generation)

    # build model
    model_equation = "y = intercept  + genotypes";
    model = build_model(model_equation);

    out   = runMCMC(model,phenotypes,chain_length=chain_length,single_step_analysis=true,pedigree=pedigree);
end
