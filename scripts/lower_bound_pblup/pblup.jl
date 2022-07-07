# this file is to run JWAS
cd("/Users/tianjing/Library/CloudStorage/Box-Box/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/pblup")
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)
ntest=100
chain_length=10000
#read pedigree
pedfile    = "/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/pedigree.txt"
pedigree   = get_pedigree(pedfile,separator=",",header=true)
# Pedigree information:
# #individuals: 6473
# #sires:       1011
# #dams:        3102
# #founders:    1247

#read phenotype
phenofile  ="/Users/tianjing/Library/CloudStorage/Box-Box/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/y.csv"
phenotypes = CSV.read(phenofile,DataFrame,missingstring=".")
phenotypes[!,:ID] = string.(phenotypes[!,:ID])
phenotypes[!,:y]= convert(Vector{Union{Missing,Float64}}, phenotypes[!,:y]); #convert the data type to allow missing data
phenotypes[end-ntest+1:end, :y] .= missing;  #missing phenotype for individuals in testing dataset
testID=phenotypes[end-ntest+1:end,:ID]

# build model
model_equation = "y = intercept + ID";
model = build_model(model_equation);
set_random(model,"ID",pedigree);

out   = runMCMC(model,phenotypes,chain_length=chain_length)

#prediction accuracy
testIndex=findall(x -> x ∈ testID, results[:,:ID])
accuruacy_test  = cor(results[testIndex,:EBV],results[testIndex,:bv]) #0.414
@show accuruacy_test



####### draw trace plot for prediction accuracy #######
mcmc_ebv  = CSV.read("results/MCMC_samples_EBV_y.txt",DataFrame) #MCMC samples for estimated BV
ind_names = names(mcmc_ebv)
mcmc_ebv  = Matrix(mcmc_ebv)
mcmc_mean = zeros(size(mcmc_ebv))
for i in 1:1000
    mcmc_mean[i,:]=mean(mcmc_ebv[1:i,:],dims=1)
end

mcmc_df=DataFrame(mcmc_mean')
insertcols!(mcmc_df,1,:ID => ind_names)
resall=innerjoin(mcmc_df, phenotypes, on = :ID)
testIndex=findall(x -> x ∈ testID, resall[:,:ID])

accuracy_datai=ones(1000)*999
for i in 1:1000
    accuracy_datai[i]=cor(resall[testIndex,"x$i"],resall[testIndex,:bv])
end

#save results
open("accuracy.txt", "w") do io
    writedlm(io, accuracy_datai)
end
final_accuracy=round(accuracy_datai[end],digits=3)

#plot results
using Plots
myfig=plot(collect(chain_length/1000:chain_length/1000:chain_length),accuracy_datai,xlabel="iteration",ylabel="accuracy",title="accu=$final_accuracy")
savefig(myfig,"accuracy.png")
