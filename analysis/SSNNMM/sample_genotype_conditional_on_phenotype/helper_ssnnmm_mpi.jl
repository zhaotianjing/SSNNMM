#this file contains helper functions for SSNNMM

#copied from JWAS
#Goal:
#  - one iteration of Gibbs sampler for lambda version of MME (single-trait)
#Arguments:
#  - A: mme.mmeLhs matrix
#  - b: mme.mmeRhs vactor
#  - x: unknowns of interest
function Gibbs!(A,x,b,vare)
    for i = 1:size(x,1)
        if A[i,i] != 0.0 #issue70, zero diagonals in MME
            invlhs  = 1.0/A[i,i]
            μ       = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn(Float32)*sqrt(invlhs*vare) + μ
        end
    end
end


#modified from JWAS code
#Goal:
# - sample σ2_g from its posterior distribution, where g~N(0,A*σ2_g)
function sampleVCs(g,nInd,Ainv,df,scale)
    S    = g'*Ainv*g
    σ2_g = (S + df*scale)/rand(Chisq(nInd + df))
    return σ2_g
end


#copied from JWAS code
# Goal:
#  - sample variance from scaled and inverse Chi-square
function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end


#the function below using pointers is deprecated
function get_column(X,j)
    nrow,ncol = size(X)
    if j>ncol||j<0
        error("column number is wrong!")
    end
    indx = 1 + (j-1)*nrow
    ptr = pointer(X,indx)
    unsafe_wrap(Array,ptr,nrow) #pointer_to_array(ptr,nrow) #deprected in Julia 0.5
end

function get_column_ref(X)
    ncol = size(X)[2]
    xArray = Array{Array{Float32,1}}(undef,ncol)
    for i=1:ncol
        xArray[i] = get_column(X,i)
    end
    return xArray
end


function sample_w1!(Z_all,nMiddle,yobs_corr,weights_NN,σ2_yobs,σ2_weightsNN,nTrain)
      #preprocess data
      z_array_all    = get_column_ref(Z_all) #[view(Z_all, :, j) for j in 1:nMiddle] #get_column_ref(Z_all) #[Z_all[:,j] for j in 1:nMiddle]      #reference, not copy, vector of vector
      ztz_all        = [dot(zj,zj) for zj in z_array_all]     #vector of zj'zj
  
      #sample intercept μ
      yobs_corr[:]     = yobs_corr .+ weights_NN[1]
      weights_NN[1] = sum(yobs_corr)/nTrain + randn(Float32)*sqrt(σ2_yobs/nTrain)
      yobs_corr[:]     = yobs_corr .- weights_NN[1]
  
      #sample marker effects w1 (RR-BLUP)
      λ  = σ2_yobs/σ2_weightsNN
      for j=1:nMiddle
          oldAlpha  = weights_NN[j+1]
          zj        = z_array_all[j]
          zjtzj     = ztz_all[j]
  
          rhs       = dot(zj,yobs_corr) + zjtzj*oldAlpha
  
          lhs       = zjtzj + λ
          invLhs    = Float32(1.0/lhs)
          mean      = invLhs*rhs
          newAlpha  = mean + randn(Float32)*sqrt(invLhs*σ2_yobs)
  
          weights_NN[j+1] = newAlpha
          BLAS.axpy!(oldAlpha-newAlpha, zj, yobs_corr) #update yobs_corr due to update of w1
      end
end

#changed: yobs_corr_missingIndex, my_Z_t
function sampleZ!(yobs_corr_missingIndex,my_Z_t, 
                  missingIndex,
                  dn, 
                  my_w1_t,
                  my_σ2_e_t,
                  σ2_yobs;
                  update_ycorr=false)
      cn  = yobs_corr_missingIndex + my_Z_t[missingIndex]*my_w1_t                #cnj
      fn  = (cn*my_w1_t*my_σ2_e_t + dn*σ2_yobs)/(my_w1_t^2 * my_σ2_e_t + σ2_yobs) #fnj
      s2  = (my_σ2_e_t*σ2_yobs)/(my_w1_t^2 * my_σ2_e_t + σ2_yobs)                 #sj2
      my_Z_t[missingIndex] = fn + randn(Float32,length(fn))*sqrt(s2)
      if update_ycorr==true
        yobs_corr_missingIndex[:] = cn - my_Z_t[missingIndex]*my_w1_t; #update yobs_corr due to update of Z (NOTE: only consider Z in this rank for parallel computing)
      end
end
