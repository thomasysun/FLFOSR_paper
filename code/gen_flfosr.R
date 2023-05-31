#----------------------------------------------------------------------------
# Generate simulated longitudinal functional data
#----------------------------------------------------------------------------

library(fda)
library(mgcv)

source("code/helper_functions.R")

#############

# N: number of subjects
# Mi: number of replicates of each subject
# L: number of non-zero scalar covariates
# Tn: total number of timepoints
# K: number of basis functions in final resulting basis (with K-2 reparameterized P-spline functions)
# sig_noise: observation level variance
# sig_with: within-curve smooth variance
# sig_bet: between-subjects smooth variance
# sig_alpha: fixed effects variance

simfosr <- function(N, Mi, L, Tn, K = 5, sig_noise = .1, sig_with = .5, sig_bet = 1, sig_alpha = 2, seed = 11){
  
  if(length(Mi) == 1){
    Mi <- rep(Mi, N)
    Mi <<- Mi  
  }
  
  N <<- N

  tau <- seq(0, 1, by = 1/(Tn-1))
  
  b = create.bspline.basis(rangeval = c(1, Tn), breaks = seq(1, Tn, 10), norder = 4)
  Bmat = eval.basis(1:Tn, b)
  Kbmat <- ncol(Bmat)
  P <- psBasis(1:Kbmat)$P
  B <- cbind(1/sqrt(Tn), poly(tau, 1), eigen(Bmat %*% ((ginv(P)) %*% t(Bmat)), symmetric = T)$vectors[,1:(K-2)])
  
  sig_w <- sig_with
  
  sig_ga <- sig_bet
  
  sig_alpha <- sig_alpha
  
  set.seed(seed)
  alpha_true <- matrix(NA, nrow = K, ncol = L+1)

  
  alpha_true[,1] <- 1/K

  for(i in 2:(L+1)){
    
    a <- rmvn(1, c(rep(0, K)), (sig_alpha)*diag(K))
    alpha_true[,i] <- a
    
  }
  alphaf_true <- B%*%alpha_true
  
  alpha_true <<- alpha_true
  alphaf_true <<- alphaf_true

  
  w_true <<- matrix(rnorm(sum(Mi)*K, mean = 0, sqrt(sig_w)), nrow = sum(Mi), ncol = K)
  ga_true <- matrix(rnorm(N*K, mean = 0, sqrt(sig_ga)), nrow = N, ncol = K)
  gaf_true <<-  B%*%t(ga_true)
  noise <<- matrix(rnorm(sum(Mi)*Tn, mean = 0, sd = sqrt(sig_noise)), nrow = Tn, ncol = sum(Mi))
  
  X <<- cbind(rep(1, N), scale(matrix(rnorm(N*L), nrow = N, ncol = L)))
  
  # Z <<- makeZ(Mi)
  z <<- rep(1:N, Mi)
  # Yk <<- w_true + Z%*%(X%*%t(alpha_true) + ga_true)
  # Y <- B%*%t(w_true + Z%*%(ga_true))
  Y <<- B%*%t(w_true + rowrep(X%*%t(alpha_true), Mi)) + t(rowrep(t(gaf_true), Mi)) + noise
  # matplot(Y, type = "l")
  # Y <<- lapply(split(data.frame(t( B%*%t(w_true + rowrep(X%*%t(alpha_true) + ga_true, Mi)) + noise )),  rep(1:N, Mi)), t)
  
}


# simfosr(N = 20, Mi = 5, L = 10, Tn = 144, K = 5, sig_noise = .1, sig_with = .5, sig_bet = 1, sig_alpha = 2)
# matplot(Y, type = "l")
# matplot(B%*%t(rowrep(X%*%t(alpha_true), Mi)), type = "l")

