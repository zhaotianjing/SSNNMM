# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

rep          = parse(Int, ARGS[1])  #rep1,...,rep10
chain_length = 20_000
ntest        = 100 #number of testing individuals

@show rep,chain_length,ntest

mainpath="/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/all_omics/"
cd(mainpath)

################################################
# single-step NN-MM
################################################
############ READ DATA ##########
#read phenotype
phenofile="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/all_no_missing/y.rep$rep.csv"
phenotypes = CSV.read(phenofile,DataFrame);

#read middle nodes (i.e.,snps)
nodefile  = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/gene_n3534_p100.csv"
nodes     = CSV.read(nodefile,DataFrame);
p         = size(nodes,2)-1
n         = size(nodes,1)
ntrain    = n-ntest
@show ntrain
@show size(nodes),p
insertcols!(nodes,2,:y => phenotypes[:,:y]); #insert phenotype into 2nd column of nodes
for i in 2:size(nodes,2) #skip ID
    nodes[!,i]= convert(Vector{Union{Missing,Float64}}, nodes[!,i]); #convert the data type to allow missing data
end
nodes[end-ntest+1:end, :y] .= missing  # missing phenotype for individuals in testing dataset

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

out   = runMCMC(model,nodes,chain_length=chain_length,printout_model_info=false,estimate_variance=true,output_folder="nnmm.rep$rep");




################################################
# GBLUP
################################################
############ READ DATA ##########
#read partial genotype
genofile   ="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/gene_n3534_p100.csv"
genotypes  = get_genotypes(genofile,separator=',',method="RR-BLUP",quality_control=false);

#read phenotype
phenofile="/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/all_no_missing/y.rep$rep.csv"
phenotypes = CSV.read(phenofile,DataFrame);
phenotypes[!,:ID] = string.(phenotypes[!,:ID]);

phenotypes[!,:y]= convert(Vector{Union{Missing,Float64}}, phenotypes[!,:y]); #convert the data type to allow missing data
phenotypes[end-ntest+1:end, :y] .= missing;  #missing phenotype for individuals in testing dataset

# build model
model_equation = "y = intercept + genotypes";
model = build_model(model_equation);

out   = runMCMC(model,phenotypes,chain_length=chain_length,output_folder="gblup.rep$rep");

