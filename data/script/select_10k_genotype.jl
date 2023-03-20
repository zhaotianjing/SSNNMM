using CSV,DataFrames,StatsBase,DelimitedFiles

## read all genotypes
genotype_all = CSV.read("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/geno_n3534_p50436.csv",DataFrame)
@show genotype_all[1:5,1:5]


id_geno = genotype_all[:,1] #individual IDs
genotype_all_matrix=Matrix(genotype_all[:,2:end])
@show size(genotype_all_matrix)


n,p=size(genotype_all_matrix)


snp_pos=sample(1:p, 10_000, replace=false, ordered=true)
writedlm("10k_snp_pos.txt",snp_pos,',')


genotype_new_matrix=genotype_all_matrix[:,snp_pos]
@show size(genotype_new_matrix)
writedlm("geno_n3534_p10000.txt",genotype_new_matrix,',')


genotype_new_df = DataFrame(genotype_new_matrix,["m$i" for i in 1:10_000])
insertcols!(genotype_new_df, 1, :ID => id_geno)
@show genotype_new_df[1:5,1:5]
CSV.write("geno_n3534_p10000.csv",genotype_new_df)
