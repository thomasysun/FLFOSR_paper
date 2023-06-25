#----------------------------------------------------------------------------
# Fast longitudinal function-on-scalar regression (FLFOSR)
#----------------------------------------------------------------------------

library(fda)
library(mgcv)
library(spikeSlabGAM)

source("code/helper_functions.R")

### Y: Tn x M Matrix of functional observations,
# Each column corresponds to one curve from a subject.
# Each row is a measurement at a single timepoint on a common grid.
### X: M x L+1 or N x L+1 design matrix of fixed effects
# Should include column of 1s for the intercept
### z: Length M vector specifying group memberships of each curve
# e.g. if there are 20 subjects with 5 repeated measurements each, z=rep(1:20, each = 5)
### k: number of basis functions
### S: number of total MCMC iterations
### S_burn: burn first # of MCMC iterations for warm-up
### a_a, b_a, a_g, b_g, a_w, b_w: gamma hyperparameters for variance priors of alpha, gamma and omega

flfosr1 <- function(Y, X, z, k = 10, S = 2000, S_burn = S/2,
                       a_a = .1, b_a = .1, a_g = .1, b_g = .1, a_w = .1, b_w = .1){
  

  Mi <- table(z)
  MM <- sum(Mi)

  N <- length(Mi)
  
  Tn <- nrow(Y)
  
  
  # tau is the obs points (length Tn)
  tau <- seq(0, 1, by = 1/(Tn-1))
  Bmat = psBasis(1:Tn, K = k)$X
  Kbmat <- ncol(Bmat)
  P <- psBasis(1:Tn, K = k)$P

  B = cbind(1/sqrt(Tn), poly(tau, 1),
            sm(tau, K = k, rankZ = .99999999,  spline.degree = 3, diff.ord = 2, centerBase = T))
  B = B/sqrt(sum(diag(crossprod(B))))
  K <- ncol(B)
  L <- ifelse(is.matrix(X) == T, ncol(X) - 1, 0)
  Dk <- diag(crossprod(B))
  D <-  rep(Dk, rep.int(MM, ncol(B)))
  Bk <- B%*%diag(1/sqrt(Dk))
  Yk <- crossprod(Y, Bk)
  group1 <- z
  
  Nx <- nrow(X)

  if(Nx == N){
    ZX <- rowrep(X, Mi)
  }else{
    ZX <- X
  }

  ## MCMC
  
  a_alph <- a_a
  b_alph <- b_a
  a_omega <- a_w
  b_omega <- b_w
  a_ga <- a_g
  b_ga <- b_g
  
  e <- 1
  alpha <- matrix(0, nrow = L+1 , K)
  ga <- (rowsum(t(Y - Bk%*%t(ZX%*%alpha)), group1)/c(Mi))%*%Bk
  w <- t(Y - Bk%*%t(ZX%*%alpha + rowrep(ga,Mi)))%*%Bk
  
  sig_e <- .1
  sig_alpha <- c(matrix(.0001, nrow = L+1))
  sig_ga <- rowSums(ga^2)/K
  sig_w <- rowSums(w^2)/K
  
  S <- S
  
  w_post <- list()
  ga_post <- list()
  alpha_post <- list()
  sig_e_post <- rep(NA, S)
  sig_alpha_post <- matrix(NA, nrow = S, ncol = L+1)
  sig_ga_post <- matrix(NA, nrow = S, ncol = N)
  sig_w_post <- matrix(NA, nrow = S, ncol = MM)
  
  nog <- which(Mi == 1)
  
  progress <- floor(seq(1, S, length.out= 11))
  
  for(s in 1:S){

      eG <- matrix(rep(1/(sig_e + sig_w), K), ncol = K)

      
      sumeG <- rowsum(eG, group1)
      
      eGh <- eG - (rowrep(sumeG/((1/sig_ga)+sumeG), Mi) )*eG
      
      if(Nx == N){
      Q_alpha <- matrix(rep(c(diag(1/sig_alpha) + crossprod(X*((sumeG - (sumeG^2)/(1/sig_ga+sumeG))[,1]), X)), K), ncol = K)
      }else{
      Q_alpha <- matrix(rep(c(diag(1/sig_alpha) + crossprod(ZX, ZX*eGh[,1])), K), ncol = K)
      }

      l_alpha <- crossprod((ZX), Yk*eGh)
      
      alpha <- sapply(1:K, function(x) {

        if(L <= Nx){
        
          ch_Q <- chol(matrix(Q_alpha[,x], nrow=L+1, ncol=L+1))
          alpha <- backsolve(ch_Q,
                             forwardsolve(t(ch_Q), l_alpha[,x]) +
                               rnorm((L+1)))

          alpha
          
        }else{
          H <- Map('*', 1/((1/sig_ga)+sumeG[,x]), tapply(eG[,x],  group1, function(v)
            if(length(v) > 1) {outer(v, v)}
            else{v}))
          
          u <- rnorm(L+1, 0, sqrt(sig_alpha))
          delta <- rnorm(N)

          Phi <- X*(sqrt(sumeG[,x] -  sapply(H, sum)))
          v = Phi%*%u + delta

          pw <- solve(tcrossprod(Phi*rep(sqrt(sig_alpha), each = N)) + diag(N),
                      (rowsum(Yk[,x],group1)/c(Mi))*(sqrt(sumeG[,x] -  sapply(H, sum))) - v)
          alpha <- u + sig_alpha*t(Phi)%*%pw
          

          alpha
        }
      })

      Q_gaij <- rep(1/sig_ga, times = K) + rowsum(eG, group1)

      l_gaij <- rowsum(eG*(Yk - ZX%*%alpha), group1)
      ga <- matrix(rnorm(N*K, (1/Q_gaij)*l_gaij, sqrt(1/Q_gaij)), nrow = N , K)
      ga[nog,] <- 0
      

      Q_wij <- 1/sig_w + 1/sig_e

        if(Nx == N){
        w <- matrix(rnorm(MM*K, (1/(sig_e))*(Yk -   rowrep(X%*%alpha + ga, Mi))*c(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
        }else{
          w <- matrix(rnorm(MM*K, (1/(sig_e))*(Yk -   ZX%*%alpha + rowrep(ga, Mi))*c(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
        }

      sig_w <- rep(1/rgamma(N, a_omega + Mi*K/2, rate = b_omega + rowSums(rowsum(w^2, group1))/2), Mi)

      sig_ga <- rep(1/rgamma(1, a_ga + N*K/2, rate = b_ga + sum(((ga))^2)/2), N)

      sig_alpha <- 1/rgamma((L+1), a_alph + K/2, rate = b_alph + rowSums((alpha)^2)/2)
      
      if(MM > 50000){
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(B, (w + rowrep((X%*%alpha + ga), Mi)) ))^2)/2)
      }else{
        if(Nx == N){
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(Bk, (w +  rowrep(X%*%alpha + ga, Mi)) ))^2)/2)
        }else{
          sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(Bk, (w +  ZX%*%alpha + rowrep(ga, Mi)) ))^2)/2)
        }
      }

      w_post[[s]] <- w
      ga_post[[s]] <- ga
      alpha_post[[s]] <- alpha
      sig_e_post[s] <- sig_e
      sig_alpha_post[s,] <- sig_alpha
      sig_ga_post[s,] <- sig_ga
      sig_w_post[s,] <- sig_w

      if(s %in% progress){
        print(paste0("MCMC draws: [", s, "/", S, "]"))
      }
   
  }
  
  #store MCMC draws of fixed effects functions
  alphaf_post <- list()
  for(i in 1:(L+1)){
    alphaf_post[[i]] <- (Bk)%*%t(do.call(rbind, lapply(alpha_post, function(x) x[i,]))[(S_burn + 1):S,])
  }
  
  m1 <- list(X = X,
             B = Bk,
             w_post = w_post,
             ga_post = ga_post,
             alpha_post = alpha_post,
             sig_e_post = sig_e_post,
             sig_alpha_post = sig_alpha_post,
             sig_ga_post = sig_ga_post,
             sig_w_post = sig_w_post,
             alphaf_post = alphaf_post)
  
  return(m1)
  
}
