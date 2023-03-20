cd("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p10000/")
using DataFrames,CSV,Random,DelimitedFiles,Statistics,Distributions

geno_df = CSV.read("geno_n3534_p10000.csv",DataFrame)
@show geno_df[1:5,1:5]
#  Row │ ID     m1       m2       m3       m4
#      │ Int64  Float64  Float64  Float64  Float64
# ─────┼───────────────────────────────────────────
#    1 │   584      1.0      1.0      1.0      2.0
#    2 │   585      0.0      2.0      1.0      2.0
#    3 │   587      2.0      2.0      0.0      1.0
#    4 │   588      1.0      1.0      1.0      0.0
#    5 │   589      0.0      1.0      2.0      2.0

@show size(geno_df) #n=3534,p=10000+1; 1st column is ID.
ind_ID = geno_df[:,1]

geno=Matrix(geno_df[:,2:end]) #n=3534,p=50436
n,p = size(geno)
@show size(geno)

#y (h2=0.7, 0.5% QTL)
n=size(geno,1)
p=size(geno,2)
nQTL=Int(floor(p*0.005))
@show nQTL

QTL_pos=sample(MersenneTwister(1),collect(1:p),nQTL;replace=false,ordered=true) #QTL positions
QTL_effect=randn(MersenneTwister(2),nQTL) #QTL effects
writedlm("QTL_pos.txt",    QTL_pos,    ',')
writedlm("QTL_effect.txt", QTL_effect, ',')

bv=geno[:,QTL_pos] * QTL_effect #breeding values
bv=bv/std(bv)*sqrt(0.7)
@show var(bv)

y=bv+randn(MersenneTwister(3),n)*sqrt(0.3) #phenotypes
@show var(y) #~1

@show var(bv)/var(y)

y_df=DataFrame(ID=ind_ID,y=y,bv=bv)
CSV.write("y.csv",y_df)



######prepare data for replications######
#in each replicate, some random individuals in the training dataset are non-genotyped
y=CSV.read("y.csv",DataFrame)
nALL=size(y,1)
Random.seed!(1)
pctAll=[0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99] #%individuals in the training dataset without genotype data
ntest=100
ntrain=nALL-ntest
@show ntrain,ntest,nALL

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
geno_all=CSV.read("geno_n3534_p10000.csv",DataFrame)
pctAll=[0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99] #%individuals in the training dataset without genotype data
n=size(geno_all,1)
for pct in pctAll
    for rep in 1:nRep
        geno_tmp          = deepcopy(geno_all)
        missing_index     = Int.(vec(readdlm("missingindex/missingindex.pct$pct.rep$rep.txt")))
        non_missing_index = setdiff(collect(1:n),missing_index)
        geno_tmp          = geno_tmp[non_missing_index,:] #genotypes for genotyped individuals
        @show size(geno_tmp)
        @show size(geno_tmp,1)+length(missing_index)
        CSV.write("single_step_partial_genotype/partialgenotype.pct$pct.rep$rep.csv",geno_tmp)
    end
end
