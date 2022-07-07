cd("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker")
using DataFrames,CSV,Random,DelimitedFiles,Statistics,Distributions

genoALL=CSV.read("/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/geno_n3534_p50436.csv",DataFrame)
geno=genoALL[:,1:101] #use the first 1000 SNPs

ind_nameALL=geno[:,:ID]

ind_names=ind_nameALL

geno=Matrix(geno[:,2:end])
snp_names=["m$i" for i in 1:100]

geno_df=DataFrame(geno,snp_names)
insertcols!(geno_df, 1, :ID => ind_names)
geno_df
CSV.write("gene_n3534_p100.csv",geno_df) #save the genotypes for 100 SNPs

#y (h2=0.7, 5% QTL)
n=size(geno,1)
p=size(geno,2)
nQTL=5

QTL_pos=sample(MersenneTwister(1),collect(1:p),nQTL;replace=false,ordered=true) #QTL positions
QTL_effect=randn(MersenneTwister(2),nQTL) #QTL effects

bv=geno[:,QTL_pos] * QTL_effect #breeding values
bv=bv/std(bv)*sqrt(0.7)
@show var(bv)

y=bv+randn(MersenneTwister(3),n)*sqrt(0.3) #phenotypes
@show var(y) #~1

@show var(bv)/var(y)

y_df=DataFrame(ID=ind_names,y=y,bv=bv)
CSV.write("y.csv",y_df)



######prepare data for replications######
#in each replicate, some random individuals in the training dataset are non-genotyped
y=CSV.read("y.csv",DataFrame)
nALL=size(y,1)
Random.seed!(1)
pctAll=[0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99] #%individuals in the training dataset without genotype data
ntest=100

ntrain=nALL-ntest
nRep=20
for pct in pctAll
    nInd_nogeno=Int(floor(ntrain*pct))
    for rep in 1:nRep
        index_rep=sample(collect(1:ntrain),nInd_nogeno;replace=false,ordered=true) #the index for non-genotyped individuals
        @show length(index_rep)
        writedlm("missingindex/missingindex.pct$pct.rep$rep.txt",index_rep)
    end
end


#save partial genotype for conventional single-step
Random.seed!(1)
geno_all=CSV.read("gene_n3534_p100.csv",DataFrame)
pctAll=[0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99] #%individuals in the training dataset without genotype data
n=size(geno_all,1)
for pct in pctAll
    for rep in 1:nRep
        geno_tmp          = copy(geno_all)
        missing_index     = Int.(vec(readdlm("missingindex/missingindex.pct$pct.rep$rep.txt")))
        non_missing_index = setdiff(collect(1:n),missing_index)
        geno_tmp          = geno_tmp[non_missing_index,:] #genotypes for genotyped individuals
        @show size(geno_tmp)
        @show size(geno_tmp,1)+length(missing_index)
        CSV.write("singlte_step_partial_genotype/partialgenotype.pct$pct.rep$rep.csv",geno_tmp)
    end
end

