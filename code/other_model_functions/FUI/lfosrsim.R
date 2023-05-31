#' Simulate longitudinal functional response Y (matrix) and a scalar predictor X (vector)
#' 
#' @param family distribution of longitudinal functional data, including "gaussian", "binomial", "poisson".
#' @param I number of subjects.
#' @param J mean number of observations per subject.
#' @param L number of grid points on the functional domain.
#' @param beta_true true fixed effects functions.
#' @param psi_true true orthonormal functions to generate subject-specific random effects.
#' @param psi2_true true orthonormal functions to generate subject/visit-specific random effects.
#' @param SNR_B relative importance of random effects
#' @param SNR_sigma signal-to-noise ratio in Gaussian data generation
#' @export
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

lfosrsim <- function(family = "gaussian", I = 100, J = 10, L = 100, 
                      beta_true, psi_true, psi2_true = NULL, 
                      SNR_B = 1, SNR_sigma = 1){
  library(mvtnorm)
  
  ## generate number of visits for each subject from poisson distribution
  J_subj <- pmax(rpois(I, J), 1)
  
  ## generate fixed effects
  n <- sum(J_subj)
  X_des = cbind(1, rnorm(n, 0, 2))
  fixef <- X_des %*% beta_true
  
  ## generate random effects
  subj <- as.factor(rep(1:I, J_subj))
  Z_des <- model.matrix( ~ 0 + subj)
  c_true <- rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function
  b_true <- c_true %*% psi_true
  ranef = Z_des %*% b_true
  if(!is.null(psi2_true)){ ## by default do not add subject-visit random deviation
    c2_true <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
    ranef <- ranef + c2_true %*% psi2_true
  }
  
  ## generate linear predictors
  ranef <- sd(fixef)/sd(ranef)/SNR_B*ranef ## adjust for relative importance of random effects
  eta_true <- fixef + ranef
  
  ## generate longitudinal functional data
  Y_obs <- matrix(NA, n, L) 
  p_true <- plogis(eta_true)
  lam_true <- exp(eta_true)
  sd_signal <- sd(eta_true)
  for(i in 1:n){
    for(j in 1:L){
      if(family == "gaussian"){
        Y_obs[i, j] <- rnorm(1, mean = eta_true[i, j], sd = sd_signal/SNR_sigma)
      }else if(family == "binomial"){
        Y_obs[i, j] <- rbinom(1, 1, p_true[i, j])
      }else if(family == "poisson"){
        Y_obs[i, j] <- rpois(1, lam_true[i, j])
      }
    }
  }
  
  ## combine simulated data
  visit <- rep(1:I, J_subj)
  for(i in 1:I){
    visit[which(subj == i)] <- 1:J_subj[i]
  }
  dat.sim <- data.frame(ID = subj, visit = visit, X = X_des[,2], Y = I(Y_obs), eta = I(eta_true))
  
  return(dat.sim)
}
