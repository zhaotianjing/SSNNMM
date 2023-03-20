#load packages
using JWAS, LinearAlgebra, SparseArrays, Distributions, Random, DelimitedFiles, CSV, DataFrames, JLD2
Random.seed!(123)

cd("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/simulate_n3534_p10000/")
datapath="/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/simulate_n3534_p10000/"


# ################################################################################
# # Step1. calculate Ainv
# ################################################################################
include("data_prep_helper.jl")
## read genotypes
genotype = CSV.read("geno_n3534_p10000.csv",DataFrame)
id_geno  = genotype[:,:ID]                   #individual IDs
id_geno  = string.(id_geno)
@show id_geno[1:5]
nInd = length(id_geno)

## read pedigrees
pedigree    = get_pedigree("pedigree.txt",separator=",", header=true)
Ainv        = AInverse(pedigree)        #calculate inverse of A (numerator relationship matrix), 6473×6473, 34863 stored entries
Ainv = Float32.(Ainv)
id_pedigree = pedigree.IDs              #individuals in Ainv, 6473x1
@show id_pedigree[1:5]

JLD2.save("Ainv_n6473.jld2", "Ainv", Ainv)  #save


