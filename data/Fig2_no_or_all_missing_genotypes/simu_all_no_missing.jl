cd("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker")
using DataFrames,CSV,Random,DelimitedFiles,Statistics,Distributions

nRep=20 #replications
geno_df=CSV.read("gene_n3534_p100.csv",DataFrame)
geno_matrix=Matrix(geno_df[:,2:end])
ind_names=geno_df[:,1]

#y (h2=0.7, 5% QTL)
n=size(geno_matrix,1)
p=size(geno_matrix,2)
nQTL=5

#simulate phenotype for each replicate
Random.seed!(1)
for rep in 1:nRep
    QTL_pos=sample(collect(1:p),nQTL;replace=false,ordered=true)
    QTL_effect=randn(nQTL)

    bv=geno_matrix[:,QTL_pos] * QTL_effect
    bv=bv/std(bv)*sqrt(0.7)

    y=bv+randn(n)*sqrt(0.3)

    y_df=DataFrame(ID=ind_names,y=y,bv=bv)
    CSV.write("all_no_missing/y.rep$rep.csv",y_df)
end
