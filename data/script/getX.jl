using JWAS, LinearAlgebra, SparseArrays, Distributions, Random, DelimitedFiles, CSV, DataFrames,JLD2
Random.seed!(123)

cd("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p10000/")
include("data_prep_helper.jl")
## read pedigrees
pedigree    = get_pedigree("pedigree.txt",separator=",", header=true)
Ainv        = AInverse(pedigree)        #calculate inverse of A (numerator relationship matrix), 6473×6473, 34863 stored entries
Ainv = Float32.(Ainv)
id_pedigree = parse.(Int, pedigree.IDs)              #individuals in Ainv, 6473x1
nInd_pedi=length(id_pedigree)
writedlm("id_pedigree.txt",id_pedigree)
JLD2.save("Ainv_n6473.jld2", "Ainv", Ainv)  #save

# minimum(Ainv) #-1.19, NOTE: Ainv will have negative values

Random.seed!(1)
geno_all=CSV.read("geno_n3534_p10000.csv",DataFrame)
genoID=geno_all[!,1]
n=size(geno_all,1)
geno_all=nothing
pctAll=[0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99] #%individuals in the training dataset without genotype data
nRep=20

for pct in pctAll
    for rep in 1:nRep
      missing_index     = Int.(vec(readdlm("missingindex/missingindex.pct$pct.rep$rep.txt")))
      non_missing_index = setdiff(collect(1:n),missing_index)
      genotyped_ID = genoID[non_missing_index] #ID for genotyped inds

      # Create an empty incidence matrix
      X = zeros(Int32,length(genotyped_ID), nInd_pedi)

      # Loop through the genoID vector and pediID vector to fill the incidence matrix
      for i in 1:length(genotyped_ID)
            for j in 1:length(id_pedigree)
                  if genotyped_ID[i] == id_pedigree[j]
                        X[i, j] = 1
                  end
            end
      end
      X = sparse(X)

      JLD2.save("X/X.pct$pct.rep$rep.jld2", "X",X)  #save
    end
end