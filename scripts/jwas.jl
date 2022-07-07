# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

method            = ARGS[1]  # ("nnmm_pblup" "ss")
rep               = parse(Int, ARGS[2])  #rep1,...,rep10
pct               = parse(Float64, ARGS[3])
chain_length      = 20000
estimate_variance = true
ntest             = 100 #number of testing individuals

@show method,rep,pct,chain_length,ntest

mainpath="/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/"
cd(mainpath*"$pct/$method/rep$rep/")

################################################
# Single-step NN-MM
################################################
if method=="nnmm_pblup"
    ############ READ DATA ##########
    #read missing index (some individuals have no middle nodes)
    missingIndex=Int.(vec(readdlm("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/missingindex/missingindex.pct$pct.rep$rep.txt")))

    #read phenotype
    phenofile="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/y.csv"
    phenotypes = CSV.read(phenofile,DataFrame);

    #read middle nodes (i.e.,snps)
    nodefile  = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/gene_n3534_p100.csv"
    nodes  = CSV.read(nodefile,DataFrame);
    p      = size(nodes,2)-1
    @show size(nodes),p
    insertcols!(nodes,2,:y => phenotypes[:,:y]); #insert phenotype into 2nd column of nodes
    for i in 2:size(nodes,2) #skip ID
        nodes[!,i]= convert(Vector{Union{Missing,Float64}}, nodes[!,i]); #convert the data type to allow missing data
    end
    nodes[end-ntest+1:end, :y] .= missing  # missing phenotype for individuals in testing dataset
    nodes[missingIndex, 3:end] .= missing; # some individuals have no SNPs
    @show [count(ismissing,col) for col in eachcol(nodes)]

    #read pedigree
    pedfile    = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/pedigree.txt"
    pedigree   = get_pedigree(pedfile,separator=",",header=true);

    # build model
    model_equation = "y = intercept + ID";
    model = build_model(model_equation;
                        num_hidden_nodes=p,
                        latent_traits=["m$i" for i in 1:p],
                        nonlinear_function="linear",
                        is_ssnnmm=true,
                        save_middle_nodes=true,
                        save_σ2_yobs=true,
                        save_σ2_weightsNN=true);
    set_random(model,"ID",pedigree);

    out   = runMCMC(model,nodes,chain_length=chain_length,printout_model_info=false,estimate_variance=estimate_variance);
end

################################################
# conventional single-step
################################################
if method=="ss"
    ############ READ DATA ##########
    #read partial genotype
    genofile   ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/single_step_partial_genotype/partialgenotype.pct$pct.rep$rep.csv"
    genotypes  = get_genotypes(genofile,separator=',',method="RR-BLUP",quality_control=false);

    #read pedigree
    pedfile    = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/pedigree.txt"
    pedigree   = get_pedigree(pedfile,separator=",",header=true);

    #read phenotype
    phenofile  ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/y.csv"
    phenotypes = CSV.read(phenofile,DataFrame);
    phenotypes[!,:ID] = string.(phenotypes[!,:ID]);

    phenotypes[!,:y]= convert(Vector{Union{Missing,Float64}}, phenotypes[!,:y]); #convert the data type to allow missing data
    phenotypes[end-ntest+1:end, :y] .= missing;  #missing phenotype for individuals in testing dataset 

    # build model
    model_equation = "y = intercept  + genotypes";
    model = build_model(model_equation);

    out   = runMCMC(model,phenotypes,chain_length=chain_length,single_step_analysis=true,pedigree=pedigree);
end

################################################
# conventional GBLUP (genotyped individuals only)
################################################
if method=="gblup"
    ############ READ DATA ##########
    #read partial genotype
    genofile   ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/single_step_partial_genotype/partialgenotype.pct$pct.rep$rep.csv"
    genotypes  = get_genotypes(genofile,separator=',',method="RR-BLUP",quality_control=false);

    #read phenotype
    phenofile  ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/y.csv"
    phenotypes = CSV.read(phenofile,DataFrame);
    phenotypes[!,:ID] = string.(phenotypes[!,:ID]);

    phenotypes[!,:y]= convert(Vector{Union{Missing,Float64}}, phenotypes[!,:y]); #convert the data type to allow missing data
    phenotypes[end-ntest+1:end, :y] .= missing;  #missing phenotype for individuals in testing dataset

    # build model
    model_equation = "y = intercept + genotypes";
    model = build_model(model_equation);

    out   = runMCMC(model,phenotypes,chain_length=chain_length);
end
