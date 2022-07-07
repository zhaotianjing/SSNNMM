ENV["GKSwstype"]="100" #to enable saving figure in server

#this file is to calculte prediction accuracy
using DataFrames,CSV,Random,DelimitedFiles,Statistics
using Plots

Random.seed!(123)
method="ss"
main_path="/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p100_qtl5_in_marker/"
pct = parse(Float64, ARGS[1]) #e.g., 0.3
nRep=10
chain_length=20_000
ntest=100
#read true bv
bv=CSV.read("/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p100_qtl5_in_marker/y.csv",DataFrame);
bv[!,:ID] = string.(bv[!,:ID]);
@show pct
testID=bv[end-ntest+1:end,:ID]

cd(main_path*"$pct/$method")
all_accuracy=ones(nRep,2)*999

for rep in 1:nRep
    @show rep
    #read results for repi
    if ispath("rep$rep/results/EBV_y.txt")
        mcmc_ebv  = CSV.read("rep$rep/results/MCMC_samples_EBV_y.txt",DataFrame)
        ind_names = names(mcmc_ebv)
        mcmc_ebv  = Matrix(mcmc_ebv)
        mcmc_mean = zeros(size(mcmc_ebv))
        for i in 1:1000
            mcmc_mean[i,:]=mean(mcmc_ebv[1:i,:],dims=1)
        end

        mcmc_df=DataFrame(mcmc_mean')
        insertcols!(mcmc_df,1,:ID => ind_names)
        resall=innerjoin(mcmc_df, bv, on = :ID)
        testIndex=findall(x -> x ∈ testID, resall[:,:ID])
        @show length(testIndex)

        accuracy_datai=ones(1000)*999
        for i in 1:1000
            accuracy_datai[i]=cor(resall[testIndex,"x$i"],resall[testIndex,:bv])
        end

        #save results
        open(main_path*"$pct/accuracy_$method/accuracy.$method.pct$pct.rep$rep.txt", "w") do io
            writedlm(io, accuracy_datai)
        end
        final_accuracy=round(accuracy_datai[end],digits=3)

        #plot results
        myfig=plot(collect(chain_length/1000:chain_length/1000:chain_length),accuracy_datai,xlabel="iteration",ylabel="accuracy",title="pct$pct.rep$rep.accu=$final_accuracy")
        savefig(myfig,main_path*"$pct/accuracy_$method/accuracy.$method.pct$pct.rep$rep.png")

        all_accuracy[rep,:] = [rep, accuracy_datai[end]]
    end
end

writedlm(main_path*"$pct/accuracy_$method/ALL_accuracy.$method.pct$pct.txt",all_accuracy)

