#----------------------------------------------------------------------------
# Replicate the simulation results in the paper
#----------------------------------------------------------------------------

library(ggplot2)
library(ggh4x)
library(tidyr)
library(posterior)
library(refund)

source("code/flfosr1.R")
source("code/gen_flfosr.R")
source("code/helper_functions.R")
source("code/other_model_functions/fgee.R")
source("code/other_model_functions/gibbs_mult_fpca2.R")
source("code/other_model_functions/FUI/lfosr3s.R")

# simcomp: generates custom simulated data and runs FLFOSR + competing models, then calculates RMSE/MCIW/ECP and other quantities
## Models: 
### 1 "prop": proposed model, FLFOSR
### 2 "gkg": gibbs sampler from refund package
### 3 "gkvb": variational bayes from refund package
### 4 "freq": frequentist fixed effects inference from Li et al. (2022) https://pubmed.ncbi.nlm.nih.gov/35491388/
### 5 "FUI": Fast Univariate Inference from Cui et al. (2022) https://pubmed.ncbi.nlm.nih.gov/35712524/
simcomp <- function(N, Mis, L, Tn,  Ksim = 5, kt = 10,
                    sig_noise = .1, sig_with = .5, sig_bet = 1, sig_alpha = 2, shp = c(.1,.1,.1,.1,.1,.1),
                    models = c("prop", "gkg", "gkvb", "freq", "fui"), seed = 11, S = 2000, Sburn = 1000){
  
  output <<- matrix(NA, nrow = length(seed), ncol =5*4 + 6 + 9)
  
  for(jjj in seed){
    
    simfosr(N, Mis, L, Tn, K = Ksim, sig_noise, sig_with, sig_bet, sig_alpha, seed = jjj)
    
    Yall1 <- makeYall()
    Xvars <- colnames(Yall1)[grep(pattern = "V", x = colnames(Yall1))][-1]
    
    if("prop" %in% models){
      s1t1 <- Sys.time()
      m1 <- flfosr1(Y, X, z, k = kt, S  = S, S_burn = Sburn, a_a = shp[1], b_a = shp[2], a_g = shp[3], b_g = shp[4], a_w = shp[5], b_w = shp[6])
      s1t2 <- difftime(Sys.time(), s1t1, units = "secs")
      esses1 <- apply(do.call(rbind, m1$alphaf_post), 1, ess_basic)
      essp1 <- (S - Sburn)*(as.numeric(s1t2)/esses1)
      ess1 <- mean(essp1)
      ess1l <- quantile(essp1, .05)
      ess1u <- quantile(essp1, .95)
      
      alphaf_mean <- apply(sapply(m1$alphaf_post, fosrcoef), 2, function(x) x$alphaf_mean)
      alphaf_l <- apply(sapply(m1$alphaf_post, fosrcoef), 2, function(x) x$alphaf_pci[1,])
      alphaf_u <- apply(sapply(m1$alphaf_post, fosrcoef), 2, function(x) x$alphaf_pci[2,])
      
      fmse1 <- sqrt(fmse(alphaf_true[,-1], apply(sapply(m1$alphaf_post, fosrcoef)[,-1], 2, function(x) x$alphaf_mean)))
      mciw1 <- mciw(alphaf_u[,-1], alphaf_l[,-1])
      ecp1 <- ecp(alphaf_u[,-1], alphaf_l[,-1], (alphaf_true[,-1]))
    }else{
      s1t2 <- NA
      ess1 <- NA
      ess1l <- NA
      ess1u <- NA
      
      fmse1 <- NA
      mciw1 <- NA
      ecp1 <- NA  
    }
    
    
    if("gkg" %in% models){
      s2t1 <- Sys.time()
      m2 <- gibbs_mult_fpca2(formula(paste(paste("Y", paste(Xvars, collapse = "+"),  sep = "~"),  paste(" + re(subject) + re(curve)"))),
                             data = Yall1, Kt = kt, N.iter = S, N.burn = Sburn)
      s2t2 <- difftime(Sys.time(), s2t1, units = "secs")
      m2alphas <- array(unlist(lapply((Sburn+1):S, function(x) m2$Theta%*%m2$BW[,,x])), dim = c(Tn, L + 1, length((Sburn+1):S)))
      esses2 <- apply(m2alphas[,-1,], 2, function(x) apply(x, 1, ess_basic))
      essp2 <- (S - Sburn)*(as.numeric(s2t2)/esses2)
      ess2 <- mean(essp2)
      ess2l <- quantile(essp2, .05)
      ess2u <- quantile(essp2, .95)
      
      fmse2 <- sqrt(fmse(alphaf_true[,-1], t(m2$beta.hat)[,-1]))
      mciw2 <- mciw(m2$beta.UB[-1,], m2$beta.LB[-1,])
      ecp2 <- ecp(m2$beta.UB[-1,], m2$beta.LB[-1,], t(alphaf_true[,-1]))
    }else{
      
      s2t2 <- NA
      ess2 <- NA
      ess2l <- NA
      ess2u <- NA
      
      fmse2 <- NA
      mciw2 <- NA
      ecp2 <- NA  
      
    }
    
    if("gkvb" %in% models){
      s3t1 <- Sys.time()
      m3 <- vb_mult_fpca(formula(paste(paste("Y", paste(Xvars, collapse = "+"),  sep = "~"),  paste(" + re(subject) + re(curve)"))),
                         data = Yall1, Kt = kt)
      s3t2 <- difftime(Sys.time(), s3t1, units = "secs")
      
      fmse3 <-sqrt(fmse(alphaf_true[,-1], t(m3$beta.hat)[,-1]))
      
      mciw3 <- mciw(m3$beta.UB[-1,], m3$beta.LB[-1,])
      ecp3 <- ecp(m3$beta.UB[-1,], m3$beta.LB[-1,], t(alphaf_true[,-1]))
    }else{
      s3t2 <- NA
      fmse3 <- NA
      mciw3 <- NA
      ecp3 <- NA  
    }
    
    if("freq" %in% models){
      ZX <- makeZ(Mi)%*%X
      s4t1 <- Sys.time()
      m4 <- fgee(formula = formula(paste("Y", paste(Xvars, collapse = "+"),  sep = "~")),
                 Y=Yall1$Y, Cov=as.data.frame(ZX)[,-1],
                 s=1:Tn, subjID=Yall1$subject, 
                 Tij=do.call(c, lapply(Mi, function(x) seq(1:x)))/max(Mi), numLongiPoints=max(Mi)*8,
                 corstr = "exchangeable")
      s4t2 <- difftime(Sys.time(), s4t1, units = "secs")
      
      fmse4 <- sqrt(fmse(alphaf_true[,-1], m4$beta[,-1]))
      
      mciw4 <- mciw(m4$beta[,-1] + m4$beta.se[,-1]*1.96, m4$beta[,-1] - m4$beta.se[,-1]*1.96)
      ecp4 <- ecp(m4$beta[,-1] + m4$beta.se[,-1]*1.96, m4$beta[,-1] - m4$beta.se[,-1]*1.96, alphaf_true[,-1])
    }else{
      s4t2 <- NA
      fmse4 <- NA
      mciw4 <- NA
      ecp4 <- NA  
    }
    
    if("fui" %in% models){    
      s5t1 <- Sys.time()
      m5 <- lfosr3s(formula = formula(paste(paste("Y", paste(Xvars, collapse = "+"),  sep = "~"),  paste(" + (1 | subject)"))), data = Yall1, family = "gaussian", 
                    var = TRUE, analytic = TRUE, parallel = FALSE, silent = TRUE)
      s5t2 <- difftime(Sys.time(), s5t1, units = "secs")      
      
      fmse5 <- sqrt(fmse(alphaf_true[,-1], t(m5$betaHat[-1,])))    
      mciw5 <- mciw(t(m5$betaHat[-1,]) + apply(m5$betaHat.var, 3, function(z) sqrt(diag(z)))[,-1]*1.96, t(m5$betaHat[-1,]) - apply(m5$betaHat.var, 3, function(z) sqrt(diag(z)))[,-1]*1.96)
      ecp5 <- ecp(t(m5$betaHat[-1,]) + apply(m5$betaHat.var, 3, function(z) sqrt(diag(z)))[,-1]*1.96, t(m5$betaHat[-1,]) - apply(m5$betaHat.var, 3, function(z) sqrt(diag(z)))[,-1]*1.96, alphaf_true[,-1])
    }else{
      s5t2 <- NA
      fmse5 <- NA
      mciw5 <- NA
      ecp5 <- NA  
    }
    
    output[which(seed == jjj),] <<- c(fmse1, fmse2, fmse3, fmse4, fmse5,
                                      mciw1, mciw2, mciw3, mciw4, mciw5,
                                      ecp1, ecp2, ecp3, ecp4, ecp5,
                                      ess1, ess1l, ess1u,
                                      ess2, ess2l, ess2u,
                                      s1t2, s2t2, s3t2, s4t2, s5t2,
                                      N, Mis, L, Tn, kt,
                                      sig_noise, sig_with, sig_bet, sig_alpha)
    
  }
  
  output <- as.data.frame(output)
  names(output) <- c("fmse1", "fmse2", "fmse3", "fmse4", "fmse5",
                     "mciw1", "mciw2", "mciw3", "mciw4", "mciw5",
                     "ecp1", "ecp2", "ecp3",  "ecp4", "ecp5",
                     "ess1", "ess1l", "ess1u",
                     "ess2", "ess2l", "ess2u",
                     "s1t2", "s2t2", "s3t2", "s4t2", "s5t2",
                     "N", "Mis", "L", "Tn", "kt",
                     "sig_noise", "sig_with", "sig_bet", "sig_alpha")
  output
  
}

# Replicate simulation results in the paper (may take a long time)

seeds <- 1:30

#### Accuracy comparison

sc1 <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc2 <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 10, sig_alpha = 1, seed = seeds)
sc3 <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 10, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc4 <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 10, sig_bet = 10, sig_alpha = 1, seed = seeds)
sc5 <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 1, sig_alpha = 10, seed = seeds)
sc6 <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
scs1 <- rbind(sc1,sc2,sc3,sc4,sc5,sc6)
#write.csv(scs1, file = "tables/dfsimresults_mse_050323.csv", row.names = FALSE)

#### Time comparison

# increasing N
sc21 <- simcomp(10, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc22 <- simcomp(20, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc23 <- simcomp(50, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc24 <- simcomp(100, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb", "fui"), seed = seeds)
sc25 <- simcomp(200, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb", "fui"), seed = seeds)

scs2 <- rbind(sc21, sc22, sc23, sc24, sc25)
#write.csv(scs2, file = "tables/dfsimresults_050323n.csv", row.names = FALSE)


# increasing Mi
sc31 <- simcomp(10, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc32 <- simcomp(10, 10, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc33 <- simcomp(10, 25, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc34 <- simcomp(10, 50, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb", "fui"), seed = seeds)
sc35 <- simcomp(10, 100, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb", "fui"), seed = seeds)
sc36 <- simcomp(10, 150, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb", "fui"), seed = seeds)

scs3 <- rbind(sc31, sc32, sc33, sc34, sc35, sc36)
#write.csv(scs3, file = "tables/dfsimresults_050323m.csv", row.names = FALSE)


# increasing L
sc41 <- simcomp(30, 5, 5, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc42 <- simcomp(30, 5, 10, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc43 <- simcomp(30, 5, 25, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds)
sc44 <- simcomp(30, 5, 33, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, models = c("prop", "freq", "gkvb", "gkg"),  seed = seeds)
sc45 <- simcomp(30, 5, 50, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, models = c("prop", "freq", "gkvb"), seed = seeds)
sc46 <- simcomp(30, 5, 100, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb"), seed = seeds)
sc47 <- simcomp(30, 5, 200, 144, kt=15, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1,  models = c("prop", "freq", "gkvb"), seed = seeds)

scs4 <- rbind(sc41, sc42, sc43, sc44, sc45, sc46, sc47)
#write.csv(scs4, file = "tables/dfsimresults_050323l.csv", row.names = FALSE)


################ Accuracy comparison
scs1$scen <- apply(scs1[,c("sig_noise", "sig_with", "sig_bet", "sig_alpha")], 1, function(x) paste(x, collapse = ", "))

dfs <- expand.grid(unique(scs1$scen), c("FLFOSR", "refund:Gibbs", "refund:VB", "Li et al. 2022", "FUI"),  c("RMSE", "MCIW", "ECP"))
dfs <- dfs[rep(seq_len(nrow(dfs)), each = length(seeds)), ]
dfs$val <- c(as.matrix(scs1[,1:(3*5)]))

e <- ggplot(dfs, aes(x = Var2, y = val))
e2 <- e + geom_boxplot(
  aes(fill = Var3),
  position = position_dodge(0.9))   + 
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1), legend.position = "none", 
        axis.title.x = element_blank(),
        text = element_text(size = 16),
        axis.title.y = element_blank()) + geom_hline(yintercept=.95, lty=1, 
                                                     color = "red", linewidth=.5) +
  ggtitle(expression(sigma[epsilon]^"2"*"*, "*sigma[omega]^"2"*"*, "*sigma[gamma]^"2"*"*, "*sigma[alpha]^"2"*"*")) 

e2  + facet_grid(Var3~Var1, scales = "free") + scale_y_facet(ROW == 1, limits = c(0, .25)) +
  scale_y_facet(ROW == 2, limits = c(0, .8)) +
  scale_y_facet(ROW == 3, limits = c(.2, 1))

############### Time comparisons
### MCMC efficiency

scs1$effr1 <- scs1$s1t2/scs1$ess1
scs1$effr2 <- scs1$s2t2/scs1$ess2


dfs <- expand.grid(unique(scs1$scen), c("FLFOSR", "refund:Gibbs"),  c("ARE"))
dfs <- dfs[rep(seq_len(nrow(dfs)), each = length(seeds)), ]
dfs$val <- c(scs1$effr1, scs1$effr2)

e <- ggplot(dfs, aes(x = Var2, y = val))
e2 <- e + geom_boxplot(
  aes(fill = Var2),
  position = position_dodge(0.9))   + 
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1), legend.position = "none", 
        axis.title.x = element_blank(),
        text = element_text(size = 16)) + geom_hline(yintercept=1, lty=1, 
                                                     color = "black", linewidth=.5) +
  labs(y = expression(bar(N)[eff]/N)) +
  geom_hline(yintercept=0, lty=1, 
             color = "black", linewidth=.5) +
  ggtitle(expression(sigma[epsilon]^"2"*"*, "*sigma[omega]^"2"*"*, "*sigma[gamma]^"2"*"*, "*sigma[alpha]^"2"*"*")) 

e2  + facet_grid(cols = vars(Var1), scales = "free") + scale_y_facet(ROW == 1, limits = c(0, 1))

### N

scs2$effr1 <- scs2$s1t2/scs2$ess1
scs2$effr2 <- scs2$s2t2/scs2$ess2

scs2 <- apply(scs2, 2, function(x) colMeans(matrix(x, nrow = 30)))

scs2 <- as.data.frame(scs2[, c("ess1","effr1", "ess2","effr2", "s3t2", "s4t2", "s5t2", 
                               "N", "Mis", "L", "Tn", "kt")])

long <- scs2 %>% 
  pivot_longer(
    cols = c("ess1", "ess2", "s3t2", "s4t2", "s5t2"), 
    names_to = "cats",
    values_to = "value"
  )
long$cats2 <- factor(long$cats)
levels(long$cats2) <- c("FLFOSR", "refund:Gibbs", "refund:VB", "Li et al. 2022", "FUI")

g <- ggplot(data = long, aes(y = value, x = N, group = cats2, color = cats2)) +
  geom_line(aes(linetype=grepl("ess", cats)), size = 2) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  geom_point(aes(shape = cats2), size = 4) +
  scale_shape_manual(values = c("FLFOSR" = 15, "refund:Gibbs" = 17, "refund:VB" = 17, "Li et al. 2022" = 19 , "FUI" = 19)) +
  scale_y_continuous(name = "Time (s)", breaks = pretty(long$value, n = 10)) +
  scale_x_continuous(name = "Number of subjects (n)") +
  theme_bw()  + guides(linetype = "none") +
  theme(legend.title = element_blank(), legend.position = c(.67, .80), legend.text = element_text(size=24),
        legend.key.size = unit(1, 'cm'),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.background = element_rect(fill="white",
                                         size=.1, linetype="solid", 
                                         colour ="black"),
        text = element_text(size = 24))
g


### m

scs3$effr1 <- scs3$s1t2/scs3$ess1
scs3$effr2 <- scs3$s2t2/scs3$ess2

scs3 <- apply(scs3, 2, function(x) colMeans(matrix(x, nrow = 30)))

scs3 <- as.data.frame(scs3[, c("ess1","effr1", "ess2","effr2", "s3t2", "s4t2", "s5t2", 
                               "N", "Mis", "L", "Tn", "kt")])

long <- scs3 %>% 
  pivot_longer(
    cols = c("ess1", "ess2", "s3t2", "s4t2", "s5t2"), 
    names_to = "cats",
    values_to = "value"
  )
long$cats2 <- factor(long$cats)
levels(long$cats2) <- c("FLFOSR", "refund:Gibbs", "refund:VB", "Li et al. 2022", "FUI")

g2 <- ggplot(data = long, aes(y = value, x = Mis, group = cats2, color = cats2)) +
  geom_line(aes(linetype=grepl("ess", cats)), size = 2) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  geom_point(aes(shape = cats2), size = 4) +
  scale_shape_manual(values = c("FLFOSR" = 15, "refund:Gibbs" = 17, "refund:VB" = 17, "Li et al. 2022" = 19 , "FUI" = 19)) +
  scale_y_continuous(name = "Time (s)", breaks = pretty(long$value, n = 10)) +
  scale_x_continuous(name = "Number of replicates per subject (m)") +
  theme_bw()  + guides(linetype = "none") +
  theme(legend.title = element_blank(), legend.position = c(.67, .80), legend.text = element_text(size=24),
        legend.key.size = unit(1, 'cm'),
        legend.spacing.y = unit(0, "mm"),
        legend.background = element_rect(fill="white",
                                         size=.1, linetype="solid", 
                                         colour ="black"),
        text = element_text(size = 24))
g2

### L

scs4$effr1 <- scs4$s1t2/scs4$ess1
scs4$effr2 <- scs4$s2t2/scs4$ess2

scs4 <- apply(scs4, 2, function(x) colMeans(matrix(x, nrow = 30)))

scs4 <- as.data.frame(scs4[, c("ess1","effr1", "ess2","effr2", "s3t2", "s4t2", "s5t2", 
                               "N", "Mis", "L", "Tn", "kt")])

long <- scs4 %>% 
  pivot_longer(
    cols = c("ess1", "ess2", "s3t2", "s4t2", "s5t2"), 
    names_to = "cats",
    values_to = "value"
  )
long$cats2 <- factor(long$cats)
levels(long$cats2) <- c("FLFOSR", "refund:Gibbs", "refund:VB", "Li et al. 2022", "FUI")

g3 <- ggplot(data = long, aes(y = value, x = L, group = cats2, color = cats2)) +
  geom_vline(xintercept = 25, linetype="dotted", 
             color = "black", size=1.5) +
  geom_line(aes(linetype=grepl("ess", cats)), size = 2) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  geom_point(aes(shape = cats2), size = 4) +
  scale_shape_manual(values = c("FLFOSR" = 15, "refund:Gibbs" = 17, "refund:VB" = 17, "Li et al. 2022" = 19 , "FUI" = 19)) +
  scale_y_continuous(name = "Time (s)", breaks = pretty(long$value, n = 10)) +
  scale_x_continuous(name = "Number of predictors (L)") +
  theme_bw()  + guides(linetype = "none") +
  theme(legend.title = element_blank(), legend.position = c(.57, .80), legend.text = element_text(size=24),
        legend.key.size = unit(1, 'cm'),
        legend.spacing.y = unit(0, "mm"),
        legend.background = element_rect(fill="white",
                                         size=.1, linetype="solid", 
                                         colour ="black"),
        text = element_text(size = 24))
g3


#### Sensitivity analysis

sc1h <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,5,1,5,1))
sc2h <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 10, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,5,1,5,1))
sc3h <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 10, sig_bet = 1, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,5,1,5,1))
sc4h <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 10, sig_bet = 10, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,5,1,5,1))
sc5h <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 1, sig_alpha = 10, seed = seeds, models = "prop", shp = c(.1,.1,5,1,5,1))
sc6h <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,5,1,5,1))
scs1h <- rbind(sc1h,sc2h,sc3h,sc4h,sc5h,sc6h)
# write.csv(scs1h, file = "tables/dfsimresults_mse_050323_senshigh.csv", row.names = FALSE)

sc1l <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,.005,.001,.005,.001))
sc2l <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 10, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,.005,.001,.005,.001))
sc3l <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 10, sig_bet = 1, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,.005,.001,.005,.001))
sc4l <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 10, sig_bet = 10, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,.005,.001,.005,.001))
sc5l <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 1, sig_with = 1, sig_bet = 1, sig_alpha = 10, seed = seeds, models = "prop", shp = c(.1,.1,.005,.001,.005,.001))
sc6l <- simcomp(20, 5, 5, 144, kt=12, sig_noise = 10, sig_with = 1, sig_bet = 1, sig_alpha = 1, seed = seeds, models = "prop", shp = c(.1,.1,.005,.001,.005,.001))
scs1l <- rbind(sc1l,sc2l,sc3l,sc4l,sc5l,sc6l)
# write.csv(scs1l, file = "tables/dfsimresults_mse_050323_senslow.csv", row.names = FALSE)


#########
scs <- scs1
scs[,c(1,6,11) + 1] <- scs1h[,c(1,6,11)]
scs[,c(1,6,11) + 2] <- scs1l[,c(1,6,11)]
scs <- scs[,-c(4,5,9,10,14,15)]

scs1$scen <- apply(scs1[,c("sig_noise", "sig_with", "sig_bet", "sig_alpha")], 1, function(x) paste(x, collapse = ", "))

dfs <- expand.grid(unique(scs1$scen), c("Original", "High", "Low"),  c("RMSE", "MCIW", "ECP"))
dfs <- dfs[rep(seq_len(nrow(dfs)), each = length(seeds)), ]
dfs$val <- c(as.matrix(scs1[,1:(3*3)]))

e <- ggplot(dfs, aes(x = Var2, y = val))
e2 <- e + geom_boxplot(
  aes(fill = Var3),
  position = position_dodge(0.9))   + 
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1), legend.position = "none", 
        axis.title.x = element_blank(),
        text = element_text(size = 16),
        axis.title.y = element_blank()) + geom_hline(yintercept=.95, lty=1, 
                                                     color = "red", linewidth=.5) +
  ggtitle(expression(sigma[epsilon]^"2"*"*, "*sigma[omega]^"2"*"*, "*sigma[gamma]^"2"*"*, "*sigma[alpha]^"2"*"*")) 

e2  + facet_grid(Var3~Var1, scales = "free") + scale_y_facet(ROW == 1, limits = c(0, .2)) +
  scale_y_facet(ROW == 2, limits = c(0, .8)) +
  scale_y_facet(ROW == 3, limits = c(.3, 1))








