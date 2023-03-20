# load packages
using MPI
using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using InteractiveUtils

# set running parameters in all ranks
rep         = parse(Int, ARGS[1])  #rep1,...,rep20
pct         = parse(Float64, ARGS[2])
method      = ARGS[3]

# load helper functions
include("/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/$method/helper_ssnnmm_mpi.jl") #Gibbs(),sample_znj(),sampleVCs(),sample_variance()

cd("/group/qtlchenggrp/tianjing/singlestep_nnlmm/pig_n3534_p10000/$pct/$method/rep$rep/")

@show method

nIteration  = 2000    #number of MCMC iterations
outFreq     = 2       #Int(nIteration/1000) #frequency of saving MCMC samples
seed        = 123     #user-defined seed

center_z    = true    #center genotypes (middle layer)
fix_z_var   = true    #fix variance components of genotypes (middle layer)

@show rep,pct,nIteration,outFreq,seed,center_z,fix_z_var

data_path   = "/group/qtlchenggrp/tianjing/singlestepdata/pig_cleveland/simulate_n3534_p10000/"  #path of data

function t(;rep=rep,pct=pct,nIteration=nIteration,outFreq=outFreq,seed=seed,
	        data_path=data_path,
	        center_z=center_z, fix_z_var=fix_z_var,
	        save_Z_mean=true,save_U_mean=true)
	MPI.Init()
	comm         = MPI.COMM_WORLD
	my_rank      = MPI.Comm_rank(comm) #current rank, e.g., 0/1/2/3, root=0
	cluster_size = MPI.Comm_size(comm) #number of all processes, e.g., 4
	@show my_rank,cluster_size
	MPI.Barrier(comm)

	nMiddle     = 10000                #number of middle nodes
	nInd        = 3534                 #number of all individuals (train+test)

	# set seed in different rank
	Random.seed!(seed+my_rank)
	MPI.Barrier(comm)

	############################################################################
	# read data in all ranks
	############################################################################
	y_df          = CSV.read(data_path*"y.csv",DataFrame)   #output layer: phenotypes for all individuals, vector of length n
	y             = Float32.(y_df[:,:y])
	y_df = nothing #for GC

	missingIndex  = Int32.(vec(readdlm(data_path*"missingindex/missingindex.pct$pct.rep$rep.txt")))  #index of non-genotyped individuals, assume the genotype of an individual is all known or all missing
	gIndex        = setdiff(collect(1:nInd),missingIndex)  #index of genotyped individuals in 3k inds
	numInd_g      = length(gIndex) #number of genotyped individuals
	trainIndex    = collect(1:3434) #index of training individuals in 3k inds (phenotypes of testing individuals are missing)

	Ainv          = JLD2.load(data_path*"Ainv_n6473.jld2")["Ainv"]     #Ainv for all individuals in pedigree, 6473-by-6473
	X             = JLD2.load(data_path*"X/X.pct$pct.rep$rep.jld2")["X"] #incidence matrix for PBLUP, #genotypedInd-by-6473
	ng_index_among_pedi=JLD2.load(data_path*"ng_index_among_pedi/ng_index_among_pedi.pct$pct.rep$rep.jld2")["ng_index_among_pedi"] #index for ng inds in 6k inds
	@show size(X),numInd_g

	nInd_pedi     = size(Ainv,1) #number of inds in pedigree, 6473

	my_Z          = readdlm(data_path*"cluster_size$cluster_size/my_Z.rank$my_rank.txt",',') #read small Z in current rank, 3k-by-smallp
	my_Z          = Float32.(my_Z)

	nInd          = Int32(size(my_Z,1))  #number of individuals, 3534
	my_batch_size = Int32(size(my_Z,2))  #number of middle nodes (i.e., SNPs) in current rank,e.g., 500 SNPs
	my_Z          = convert(Matrix{Union{Missing,Float32}}, my_Z) #set Z as missing for non-genotyped individuals
	my_Z[missingIndex,:] .= missing

	if center_z==true #center known genotypes zg
		colmean=Float32.([mean(skipmissing(my_Z[:,i])) for i in 1:size(my_Z, 2)]')  #row vector
		my_Z = my_Z .- colmean
	end

	############################################################################
	# set model parameters in all ranks
	############################################################################
	varz = Float32.([var(skipmissing(my_Z[:,i])) for i in 1:size(my_Z, 2)]) #var(zj), length of my_batch_size
	println("min,mean,max var(z):", minimum(varz),",",mean(varz),",",maximum(varz))

	my_σ2_e    = Float32.(varz*0.001)  #residual variance in PBLUP, length of my_batch_size
	my_σ2_g    = Float32.(varz*0.999)  #genetic variance in PBLUP, assume h2=0.5, length of my_batch_size
	my_h2 = my_σ2_g ./ (my_σ2_g+my_σ2_e) #snp h2
	println("min,mean,max SNP σ2_g:", round(minimum(my_σ2_g),digits=2),",",round(mean(my_σ2_g),digits=2),",",round(maximum(my_σ2_g),digits=2),
		    "; min,mean,max SNP σ2_e:", round(minimum(my_σ2_e),digits=2),",",round(mean(my_σ2_e),digits=2),",",round(maximum(my_σ2_e),digits=2),
			"; min,mean,max SNP h2:", round(minimum(my_h2),digits=2),",",round(mean(my_h2),digits=2),",",round(maximum(my_h2),digits=2))
	
	if fix_z_var==false #sample var z
		df_e       = Int32(4)                          #degrees of freedom for σ2_e in PBLUP
		df_g       = Int32(5)                          #degrees of freedom for σ2_g in PBLUP
		my_scale_e = Float32.(my_σ2_e*(df_e-2)/df_e)   #scale of residual variance in PBLUP, length of my_batch_size
		my_scale_g = Float32.(my_σ2_g*(df_g-2))        #scale of genetic variance in PBLUP, length of my_batch_size
	end

	#initialize missing values using phenotype to convert datatype
	my_Z[missingIndex,:] = repeat(y[missingIndex],my_batch_size)
	my_Z                 = convert(Matrix{Float32}, my_Z) #convert datatype to avoid error in MPI
	my_Z_array           = get_column_ref(my_Z);

	XtX        = X'X                                    #LHS without covariances part, sparse, 6k-by-6k

	my_sol     = [zeros(Float32,nInd_pedi) for t in 1:my_batch_size] #vector of vector, each element is the [μ0;g0] in a PBLUP, length of my_batch_size
    
    if save_U_mean==true
    	u_all_mean = zeros(Float32, nInd_pedi, my_batch_size) #6k-by-#my_batch_size
	end
	
	σ2_yobs    = Float32(var(y[trainIndex])*0.5)  #initialize residual variance for observed phenotypes

	batch_size_all = MPI.Gather(my_batch_size, 0, comm) #send my_batch_size to root

	############################################################################
	# initial data and buffer in rank 0
	############################################################################
	if my_rank == 0
		nOutput         = Int32(floor(nIteration/outFreq))
		outfile1        = open("ssnnmm.w1.pct$pct.rep$rep.niter$nIteration.txt"       ,"w")
		if fix_z_var==false #sample z var
			outfile2        = open("ssnnmm.varg.pct$pct.rep$rep.niter$nIteration.txt"     ,"w") #varg for marker
			outfile3        = open("ssnnmm.vare.pct$pct.rep$rep.niter$nIteration.txt"     ,"w") #vare for marker
		end
		outfile4        = open("ssnnmm.vareyobs.pct$pct.rep$rep.niter$nIteration.txt" ,"w")
		outfile5        = open("ssnnmm.varw1.pct$pct.rep$rep.niter$nIteration.txt"    ,"w")
		outfile6        = open("mcmc.u.rank$my_rank.pct$pct.rep$rep.niter$nIteration.txt"    ,"w")

		yobs            = y[trainIndex]
		nTrain          = Int32(length(trainIndex))           #number of individuals in the training dataset

		weights_NN      = [Float32(mean(yobs)); zeros(Float32,nMiddle)]         #initialize μ1,w1

	    yobs_corr       = yobs .- weights_NN[1]        #initialize residual of phenotype, length of nTrain, y-μ-Zw1, when weights_NN[2:end]=0, and weights_NN[1]=mean(yobs)

		σ2_weightsNN       = Float32(var(yobs*0.5)/nMiddle)  #variance of neural network weight between middle and output layers, assume h2=0.5
		df_σ2_weightsNN    = Int32(4.0)
		scale_σ2_weightsNN = Float32(σ2_weightsNN*(df_σ2_weightsNN-2)/df_σ2_weightsNN)

		df_σ2_yobs         = Int32(4.0)
		scale_σ2_yobs      = Float32(σ2_yobs*(df_σ2_yobs-2)/df_σ2_yobs)

		#### save results in rank0
		Z_all           = zeros(Float32, nTrain, nMiddle) #nTrain-by-nMiddle
		counts          = batch_size_all*nTrain  # number of elements in each rank
		Z_all_vbuf      = VBuffer(Z_all, counts) # VBuffer for gather
		
		if save_Z_mean==true
			Z_all_mean  = zeros(Float32, nTrain, nMiddle) #for saving the posterior of Z
		end
		
		if fix_z_var==false #sample z var
			σ2_g_all        = zeros(Float32, nMiddle)
			σ2_g_all_vbuf   = VBuffer(σ2_g_all, batch_size_all) # VBuffer for saving results

			σ2_e_all        = zeros(Float32, nMiddle)
		    σ2_e_all_vbuf   = VBuffer(σ2_e_all, batch_size_all) # VBuffer for saving results
		end
	else
		Z_all_vbuf      = VBuffer(nothing)
		
		if fix_z_var==false #sample z var
			σ2_g_all_vbuf   = VBuffer(nothing)
			σ2_e_all_vbuf   = VBuffer(nothing)
		end
	end

	############################################################################
	# print information
	############################################################################
	if my_rank == 0
		println("---------------- Summary Start --------------")
		println("rep=$rep, pct=$pct, seed=$seed, nthreads=", Threads.nthreads())
		println("The number of all individuals: $nInd. (training: $nTrain, testing: ", nInd-nTrain,").")
		println("The number of non-genotyped individuals: ", length(missingIndex))
		println("The number of all middle nodes: $nMiddle.")
		println("The number of MCMC iterations: $nIteration. (outout freqency: $outFreq)")
		println("The number of elements in Ainv: ", nnz(Ainv), ", size of Ainv: ", size(Ainv),", type of Ainv: ",typeof(Ainv))
		println("The number of elements in X: ", nnz(X), "; size of X: ", size(X),", type of X: ", typeof(X))
		t_start = now()
		t_iter_old = now()
		println("Start time: ", t_start)
		println("---------------- Summary End ----------------")
	end
	MPI.Barrier(comm)
	println(" --Rank$my_rank: $my_batch_size middle nodes (i.e.,SNPs)") #print in all ranks
	MPI.Barrier(comm)

	############################################################################
	# MCMC iterations start (run in current rank)
	############################################################################
	for iter in 1:nIteration
	    for t in 1:my_batch_size #@Threads.threads
			####################################################################
			# Step1. run some PBLUP in current rank
			####################################################################
			# sample un
			zg_t     = my_Z_array[t][gIndex] #zg is z of genotyped inds, vector of length #gnotypedInd
			my_Rhs_t = X'zg_t
			λ_t      = my_σ2_e[t]/my_σ2_g[t]
			my_Lhs_t = Ainv*λ_t+XtX
			Gibbs!(my_Lhs_t,my_sol[t],my_Rhs_t,my_σ2_e[t]) #changed: sol[t]
			
			if save_U_mean==true
				u_all_mean[:,t] = u_all_mean[:,t] + (my_sol[t]-u_all_mean[:,t])*(1/iter)
			end

			if fix_z_var==false #sample z var
				# sample σ2_g
				my_σ2_g[t] = sampleVCs(my_sol[t], nInd_pedi, Ainv, df_g, my_scale_g[t])
				
				# sample σ2_e
				μ_zg_t     = X*my_sol[t] #vector of length #genotypedInd
				zg_corr_t  = zg_t - μ_zg_t
				my_σ2_e[t] = sample_variance(zg_corr_t, numInd_g, df_e, my_scale_e[t])
			end

			####################################################################
			# Step2. sample missing middle nodes (i.e., missing genotypes)
			####################################################################
			# assume h2=1, zn=un
      		my_Z_array[t][missingIndex] = my_sol[t][ng_index_among_pedi] #ng_index_among_pedi is pos of ng among 6k inds
      		
		end
		# send z to rank0 for step3 (training individual only)
		MPI.Gatherv!(my_Z[trainIndex,:], Z_all_vbuf, 0, comm) #Z_all: nTrain*nMiddle
		MPI.Barrier(comm)
		########################################################################
	    # Step3. sample weights between middle layer and output layer in rank 0
		########################################################################
		if my_rank == 0
			if save_Z_mean==true
				Z_all_mean  += (Z_all-Z_all_mean)*(1/iter)
			end
			#update ycorr using all new Z
			yobs_corr = yobs .- weights_NN[1] - Z_all*weights_NN[2:end] #vector of length nTrain
			#sample w1
			sample_w1!(Z_all,nMiddle,yobs_corr,weights_NN,σ2_yobs,σ2_weightsNN,nTrain)

			#sample σ2_weightsNN (marker effect variance)
			σ2_weightsNN = Float32(sample_variance(weights_NN[2:end], nMiddle, df_σ2_weightsNN, scale_σ2_weightsNN))

			#sample σ2_yobs (residual variance of phenotype)
			σ2_yobs      = Float32(sample_variance(yobs_corr, nTrain, df_σ2_yobs, scale_σ2_yobs))
		end
		MPI.Barrier(comm) #new added
		
		########################################################################
		# Step4. save mcmc samples in rank 0
		########################################################################
		if iter%outFreq == 0
			if fix_z_var==false #sample var z
				MPI.Gatherv!(my_σ2_g, σ2_g_all_vbuf, 0, comm) # for saving results
				MPI.Gatherv!(my_σ2_e, σ2_e_all_vbuf, 0, comm) # for saving results
				MPI.Barrier(comm) #new added
			end
			if my_rank == 0
				t_iter_new  = now()
				t_iter_diff = (t_iter_new-t_iter_old).value/60000 #milliseconds to min
				t_iter_old  = t_iter_new
				println("Iteraton $iter. Time(min):", round(t_iter_diff,digits=2),
					    "; σ2_weightsNN:", σ2_weightsNN,
					    "; σ2_yobs:", σ2_yobs)
				if fix_z_var==false #sample var z
					h2 = σ2_g_all ./ (σ2_g_all+σ2_e_all)
					println("; min,mean,max σ2_g:", round(minimum(σ2_g_all),digits=2),",",round(mean(σ2_g_all),digits=2),",",round(maximum(σ2_g_all),digits=2),
					        "; min,mean,max σ2_e:", round(minimum(σ2_e_all),digits=2),",",round(mean(σ2_e_all),digits=2),",",round(maximum(σ2_e_all),digits=2),
					        "; min,mean,max h2:", round(minimum(h2),digits=2),",",round(mean(h2),digits=2),",",round(maximum(h2),digits=2))
				end
				writedlm(outfile1, weights_NN',  ',') #w1
				writedlm(outfile4, σ2_yobs,      ',') #vareyobs
				writedlm(outfile5, σ2_weightsNN, ',') #varw1
				writedlm(outfile6, vcat(my_sol[1:10]...)', ',') #same mcmc samples of u for first 10 markers
				if fix_z_var==false #sample var z
					writedlm(outfile2, σ2_g_all',    ',') #varg
					writedlm(outfile3, σ2_e_all',    ',') #vare
				end
			end
		end
		MPI.Barrier(comm) #new added
	end #MCMC iterations end

	############################################################################
	# Step5. show running time
	############################################################################
	if my_rank == 0
		t_end = now()
		t_diff = (t_end-t_start).value/60000 #milliseconds to min
		println("End time: ", t_end)
		println("Running Time (min): ", t_diff)
	end

	############################################################################
	# Step6. save some results
	############################################################################
	if save_Z_mean==true
		if my_rank == 0
			writedlm("ssnnmm.Z_train_mean.pct$pct.rep$rep.niter$nIteration.txt", Z_all_mean, ',')
		end
	end
	if save_U_mean==true
		writedlm("ssnnmm.u_mean.rank$my_rank.pct$pct.rep$rep.niter$nIteration.txt", u_all_mean, ',')
	end

	MPI.Finalize()
end


t()
