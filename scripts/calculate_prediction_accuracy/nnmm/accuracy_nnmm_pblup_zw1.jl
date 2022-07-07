ENV["GKSwstype"]="100" #to enable saving figure in server

#this file is to calculte prediction accuracy
using DataFrames,CSV,Random,DelimitedFiles,Statistics
using Plots

Random.seed!(123)
method="nnmm_pblup"
main_path="/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/"
pct = parse(Float64, ARGS[1]) #e.g., 0.3
nRep=10
chain_length=20_000
ntest=100
#read true bv
bv=CSV.read("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/y.csv",DataFrame);
bv[!,:ID] = string.(bv[!,:ID])
@show pct
testID=bv[end-ntest+1:end,:ID]
bv_test=bv[end-ntest+1:end,:bv]
@show testID
cd(main_path*"$pct/$method")
all_accuracy=ones(nRep,2)*999

geno=CSV.read("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/gene_n3534_p100.csv",DataFrame)
geno_test=Matrix(geno[end-ntest+1:end,2:end]) #genotype for testing inds
geno_test_id=string.(geno[end-ntest+1:end,:ID])
@show geno_test_id==testID
@show size(geno_test)

for rep in 1:nRep
    @show rep
    #read results for repi
    if ispath("rep$rep/results/EBV_NonLinear.txt")
        mcmc_w    = readdlm("rep$rep/results/MCMC_samples_neural_networks_bias_and_weights.txt",',') #marker effects (weights between middle and output layers)
        mcmc_mean = zeros(size(mcmc_w))  #(1000,1001)
        for i in 1:1000
            mcmc_mean[i,:]=mean(mcmc_w[1:i,:],dims=1)
        end

        accuracy_datai=ones(1000)*999
        for i in 1:1000
            ebvi = geno_test*vec(mcmc_mean[i,2:end])  #Z*W1
            accuracy_datai[i]=cor(ebvi,bv_test)
        end

        #save results
        open(main_path*"$pct/accuracy_nnmm_pblup_zw1/accuracy.nnmm_pblup_zw1.pct$pct.rep$rep.txt", "w") do io
            writedlm(io, accuracy_datai)
        end
        final_accuracy=round(accuracy_datai[end],digits=3)

        #plot results
        myfig=plot(collect(chain_length/1000:chain_length/1000:chain_length),accuracy_datai,xlabel="iteration",ylabel="accuracy",title="pct$pct.rep$rep.accu=$final_accuracy")
        savefig(myfig,main_path*"$pct/accuracy_nnmm_pblup_zw1/accuracy.nnmm_pblup_zw1.pct$pct.rep$rep.png")

        all_accuracy[rep,:] = [rep, accuracy_datai[end]]
    end
end

writedlm(main_path*"$pct/accuracy_nnmm_pblup_zw1/ALL_accuracy.nnmm_pblup_zw1.pct$pct.txt",all_accuracy)

