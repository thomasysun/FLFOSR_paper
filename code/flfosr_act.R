#----------------------------------------------------------------------------
# Replicate application results in the paper
#----------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(tidyr)
library(patchwork)
library(MetBrewer)

source("code/process_nhanes_data_ts.R")
source("code/flfosr1.R")
source("code/helper_functions.R")

X[,c(2,7,10,18,19,20)] <- scale(X[,c(2,7,10,18,19,20)])

Y <- Y[,order(W)]

Y <- as.matrix(Y)

L <- dim(X)[2] - 1
Tn <- dim(Y)[1]
Mi <- table(W)
N <- length(Mi)

S <- 2000

s1t1 <- Sys.time()
m1 <- flfosr1(Y, X, z = W, k = 10, S  = S)
s1t <- difftime(Sys.time(), s1t1, units = "secs")      

### Plots the fixed effects functions
fosrcoefplots <- function(alpha_post, B, index = 1:(L+1), maintitle = ""){
  
  S <- length(alpha_post)
  L <- nrow(alpha_post[[1]]) - 1
  
  plot_list <- list()
  par(mfrow=c(2,2))
  for(i in 1:length(index)){
    
    alphaf_post <- B%*%t(do.call(rbind, lapply(alpha_post, function(x) x[index[i],]))[(S/2):S,])
    alphaf_pci <- apply(alphaf_post, 1, function(x) quantile(x, c(.025,.975)))
    alphaf_mean <- B%*%(colMeans(do.call(rbind, lapply(alpha_post, function(x) x[index[i],]))[(S/2):S,]))
    # matplot(cbind(alphaf_mean, t(alphaf_pci)), type="l", lty=c(1,2,2), col = c(1,2,2),
    #         main = colnames(X)[i])
    # #,ylim = c(-1,1))
    # 
    # abline(h=0)
    df.act.m <- data.frame(time = seq(1,nrow(alphaf_post)), 
                           mean = alphaf_mean, lower = alphaf_pci[1,], upper = alphaf_pci[2,])
    # names(df.act.m) <-c("time", "SEQN_31128", "SEQN_31193")
    df.act.m$time <- factor(df.act.m$time)
    df.act.m$time <- as.numeric(df.act.m$time)
    p.act <- ggplot(df.act.m, aes(x = time))+ 
      geom_line(aes(y = mean),size = 2, colour = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      scale_x_continuous( breaks =  seq(0, nrow(alphaf_post), len=6)[] + 1,
                          labels = c("04:00", "08:00", "12:00", "16:00", "20:00", "00:00"))+
      scale_colour_met_d("Hokusai2")+
      theme_bw()+ theme( legend.position="none",
                         axis.text.x = element_text( size = 9, colour = "black", angle = 0, vjust = 0.0, hjust=0.5),
                         axis.text.y = element_text( size = 10, colour = "black"),
                         axis.title = element_text( size = 10, colour = "black", face = "bold"))+
      xlab("Time of Day")+ylab("") + 
      geom_hline(yintercept=0) +
      if(maintitle == ""){
        ggtitle(as.name(colnames(X)[index[i]]))
      }else{
        ggtitle(as.name(maintitle[i]))
      }

    plot_list[[i]] <- p.act
    
  }
  print(ggarrange(plotlist = plot_list, nrow = 4, ncol = 4))
  par(mfrow=c(1,1))
}
titles <- c("Intercept", "BMI", "RaceBlack", "RaceHispanic", "RaceOther", "GenderFemale",
  "Age", "Education = HS", "Education > HS", "DrinksPerWeek", "Formerly Smoke:Yes", "Currently Smoke: Yes",
  "Diabetes: Yes", "CHF: Yes", "CHD: Yes", "Cancer: Yes", "Stroke: Yes", "HDL Cholesterol",
  "Total Cholesterol", "Systolic Blood Pressure", "Weekend: Yes")

fosrcoefplots(m1$alpha_post, m1$B, maintitle = titles)

fosrcoefplots(m1$alpha_post, m1$B, index = c(1:21)[-c(1,2,7,9,12,13,18,21)],
              maintitle = titles[-c(1,2,7,9,12,13,18,21)]) #9x6 in

fosrcoefplots(m1$alpha_post, m1$B, index = c(1,2,7,9,12,13,18,21), 
              maintitle = c("Intercept", "BMI", "Age", "Education > HS", "Currently Smoke: Yes",
                            "Diabetes: Yes", "HDL Cholesterol", "Weekend: Yes")) #11x5.5 in


### ESS

esses1 <- apply(do.call(rbind, m1$alphaf_post), 1, ess_basic)
essp1 <- (S - S/2)*(as.numeric(s1t)/esses1)
ess1 <- mean(essp1)
ess1l <- quantile(essp1, .05)
ess1u <- quantile(essp1, .95)

###

Yall1 <- makeYall()
Xvars <- colnames(X)[-1]

###### Prediction plots
# NOTE: need to run model without standardizing X for this to make sense
subjw_post <- array(NA, dim = c(sum(Mi), ncol(m1$B),S/2))
Xnew <- X
row.names(Xnew) <- 1:nrow(Xnew)
Xnew[,"WeekendTRUE"] <- 0 #just assume it's not the weekend
newrow <- c(1, rep(0, L))
newrow[c(2,7,10,18,19,20)] <- apply(Xnew[,c(2,7,10,18,19,20)], 2, mean) 
Xnew[1:3,] <- rbind(newrow, newrow, newrow)
Xnew[1,"Age"] <- 35
Xnew[2,"Age"] <- 45
Xnew[3,"Age"] <- 55
for(s in 1:(S/2)){
  
  subjw_post[,,s] <- Xnew%*%m1$alpha_post[[s+S/2]] 
  
}

alphaf_mean <- m1$B%*%t(apply(subjw_post[1:3,,],c(1,2),mean))
alphaf_pcil <- matrix(apply(apply(subjw_post[1:3,,], 3, function(z) m1$B%*%t(z)), 1, function(y) quantile(y, .025)), nrow=nrow(Y), ncol=3)
alphaf_pciu <- matrix(apply(apply(subjw_post[1:3,,], 3, function(z) m1$B%*%t(z)), 1, function(y) quantile(y, .975)), nrow=nrow(Y), ncol=3)

df.act.m <- data.frame(alphaf_mean)
df.act.m$time <- seq(1,nrow(alphaf_mean))
df.act.m$time <- factor(df.act.m$time)
df.act.m$time <- as.numeric(df.act.m$time)
colnames(df.act.m)[1:3] <- c("35","45","55")
df.act.m <- gather(df.act.m, Age, Activity, c("35","45","55"), factor_key = T)

df.act.l <- data.frame(alphaf_pcil)
colnames(df.act.l)[1:3] <- c("35","45","55")
df.act.l <- gather(df.act.l, Age, Activity.l, c("35","45","55"), factor_key = T)

df.act.u <- data.frame(alphaf_pciu)
colnames(df.act.u)[1:3] <- c("35","45","55")
df.act.u <- gather(df.act.u, Age, Activity.u, c("35","45","55"), factor_key = T)

df.act <- bind_cols(df.act.m, Activity.l = df.act.l[,2], Activity.u = df.act.u[,2])

p.act <- ggplot(df.act, aes(y = Activity, x = time, color = Age)) + 

  geom_ribbon(aes(ymin = Activity.l, ymax = Activity.u, fill = Age), alpha = .1, lty = 2) + 
  geom_line(size = 1.5) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  scale_x_continuous( breaks =  seq(0, nrow(alphaf_mean), len=6)[] + 1,
                      labels = c("4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM"))+
  scale_colour_met_d("Derain", direction = -1)+
  scale_fill_met_d("Derain", direction = -1)+
  theme_bw() + theme( legend.position=c(0.90, 0.90),
                     axis.text.x = element_text( size = 14, colour = "black", angle = 90, vjust = 0.5, hjust=1),
                     axis.text.y = element_text( size = 14, colour = "black"),
                     axis.title = element_text( size = 15, colour = "black", face = "bold"),
                     legend.text=element_text(size=14),
                     legend.title=element_text(size=14),
                     legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  xlab("") + 
  ggtitle("") +
  geom_hline(yintercept=0)
p.act #7x5 in
###

subjw_post <- array(NA, dim = c(sum(Mi), ncol(m1$B),S/2))
Xnew <- X
row.names(Xnew) <- 1:nrow(Xnew)
Xnew[,"WeekendTRUE"] <- 0 #just assume it's not the weekend
newrow <- c(1, rep(0, L))
newrow[c(2,7,10,18,19,20)] <- apply(Xnew[,c(2,7,10,18,19,20)], 2, mean) 
Xnew[1:3,] <- rbind(newrow, newrow, newrow)
Xnew[1,"HDL_Cholesterol"] <- 40
Xnew[2,"HDL_Cholesterol"] <- 65
Xnew[3,"HDL_Cholesterol"] <- 90
for(s in 1:(S/2)){
  
  subjw_post[,,s] <- Xnew%*%m1$alpha_post[[s+S/2]] 
  
}

alphaf_mean <- m1$B%*%t(apply(subjw_post[1:3,,],c(1,2),mean))
alphaf_pcil <- matrix(apply(apply(subjw_post[1:3,,], 3, function(z) m1$B%*%t(z)), 1, function(y) quantile(y, .025)), nrow=nrow(Y), ncol=3)
alphaf_pciu <- matrix(apply(apply(subjw_post[1:3,,], 3, function(z) m1$B%*%t(z)), 1, function(y) quantile(y, .975)), nrow=nrow(Y), ncol=3)

df.act.m <- data.frame(alphaf_mean)
df.act.m$time <- seq(1,nrow(alphaf_mean))
df.act.m$time <- factor(df.act.m$time)
df.act.m$time <- as.numeric(df.act.m$time)
colnames(df.act.m)[1:3] <- c("40","65","90")
df.act.m <- gather(df.act.m, HDL_Cholesterol, Activity, c("40","65","90"), factor_key = T)

df.act.l <- data.frame(alphaf_pcil)
colnames(df.act.l)[1:3] <- c("40","65","90")
df.act.l <- gather(df.act.l, HDL_Cholesterol, Activity.l, c("40","65","90"), factor_key = T)

df.act.u <- data.frame(alphaf_pciu)
colnames(df.act.u)[1:3] <- c("40","65","90")
df.act.u <- gather(df.act.u, HDL_Cholesterol, Activity.u, c("40","65","90"), factor_key = T)

df.act <- bind_cols(df.act.m, Activity.l = df.act.l[,2], Activity.u = df.act.u[,2])

p.act <- ggplot(df.act, aes(y = Activity, x = time, color = HDL_Cholesterol)) + 
  
  geom_ribbon(aes(ymin = Activity.l, ymax = Activity.u, fill = HDL_Cholesterol), alpha = .1, lty = 2) + 
  geom_line(size = 1.5) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  scale_x_continuous( breaks =  seq(0, nrow(alphaf_mean), len=6)[] + 1,
                      labels = c("4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM"))+
  scale_colour_met_d("Derain", direction = -1)+
  scale_fill_met_d("Derain", direction = -1)+
  labs( fill = "HDL Cholesterol \n (mg/dL)", colour = "HDL Cholesterol \n (mg/dL)", labels = "HDL Cholesterol \n (mg/dL)" ) +
  theme_bw() + theme( legend.position=c(0.87, 0.90),
                      axis.text.x = element_text( size = 14, colour = "black", angle = 90, vjust = 0.5, hjust=1),
                      axis.text.y = element_text( size = 14, colour = "black"),
                      axis.title = element_text( size = 15, colour = "black", face = "bold"),
                      legend.text=element_text(size=14),
                      legend.title=element_text(size=14),
                      legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  xlab("") + 
  ggtitle("") +
  geom_hline(yintercept=0)
p.act #7x5 in


###### Comparison with refund:Gibbs and refund:VB (refund:Gibbs takes a long time!!!)
S <- 2000
s2t1 <- Sys.time()
m2 <- gibbs_mult_fpca2(formula(paste(paste("Y", paste(Xvars, collapse = "+"),  sep = "~"),  paste(" + re(subject) + re(curve)"))),
                       data = Yall1, Kt = 10, N.iter = S, N.burn = S/2)
s2t2 <- difftime(Sys.time(), s2t1, units = "secs")

s3t1 <- Sys.time()
m3 <- vb_mult_fpca(formula(paste(paste("Y", paste(Xvars, collapse = "+"),  sep = "~"),  paste(" + re(subject) + re(curve)"))),
                   data = Yall1, Kt = 10)
s3t2 <- difftime(Sys.time(), s3t1, units = "secs")


m2alphas <- array(unlist(lapply((S/2+1):S, function(x) m2$Theta%*%m2$BW[,,x])), dim = c(120, 21, length((S/2+1):S)))
esses2 <- apply(m2alphas, c(1,2), ess_basic)
essp2 <- (S - Sburn)*(as.numeric(s2t2)/esses2)
ess2 <- mean(essp2)
ess2l <- quantile(essp2, .05)
ess2u <- quantile(essp2, .95)




###### Miscellanous plots (not in paper)
subjw_post <- array(NA, dim = c(sum(Mi),ncol(m1$B),S/2))
for(s in 1:(S/2)){
  
  subjw_post[,,s] <- X%*%m1$alpha_post[[s+S/2]] + rowrep(m1$ga_post[[s+S/2]], Mi) + m1$w_post[[s+S/2]]
  
}
subjw_pci <-
  apply(subjw_post[1:6,,], 1, function(z) apply(m1$B%*%(z), 1, function(x) quantile(x, c(.025))))
  
matplot( Y[,1:Mi[1]], type = "p")
matplot(m1$B%*%t(apply(subjw_post[1:6,,],c(1,2),mean)), add = T, type = "l", lty =1)
matplot(apply(subjw_post[1:6,,], 1, function(z) apply(m1$B%*%(z), 1, function(x) quantile(x, c(.025)))), add = T, type = "l", lty = 2)
matplot(apply(subjw_post[1:6,,], 1, function(z) apply(m1$B%*%(z), 1, function(x) quantile(x, c(.975)))), add = T, type = "l", lty = 2)


matplot( Y[,7:12], type = "p")
matplot(m1$B%*%t(apply(subjw_post[7:12,,],c(1,2),mean)), add = T, type = "l", lty =1)
matplot(apply(subjw_post[7:12,,], 1, function(z) apply(m1$B%*%(z), 1, function(x) quantile(x, c(.025)))), add = T, type = "l", lty = 2)
matplot(apply(subjw_post[7:12,,], 1, function(z) apply(m1$B%*%(z), 1, function(x) quantile(x, c(.975)))), add = T, type = "l", lty = 2)

######## EDA
Yeda <- data.frame(cbind(apply(Y, 1, function(x) mean(x)), apply(Y, 1, function(x) median(x)), apply(Y, 1, function(x) quantile(x, .75)), apply(Y, 1, function(x) quantile(x,.25))))
colnames(Yeda) <- c("Mean", "Median", "75th Quantile", "25th Quantile")
Yeda$time <- seq(1,nrow(Yeda))
Yedal <- pivot_longer(Yeda, -time, values_to = "value", names_to = "ts")
ggplot(Yedal, aes(x = time))+ 
  geom_line(aes(y = value, colour = ts),size = 1.5) +
  scale_x_continuous( breaks =  seq(0, nrow(Yeda), len=6)[] + 1,
                      labels = c("4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM"))+
  scale_colour_met_d("Derain")+
  theme_bw()+ theme( legend.position=c(0.90, 0.85),
                     legend.title=element_blank(),
                     axis.text.x = element_text( size = 13, colour = "black", angle = 90, vjust = 0.5, hjust=1),
                     axis.text.y = element_text( size = 13, colour = "black"),
                     axis.title = element_text( size = 13, colour = "black", face = "bold"),
                     legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  xlab("")+ylab("Activity") + 
  geom_hline(yintercept=0)

p1 <- ggplot(data_analysis, aes(x = BMI)) + geom_histogram()
p2 <- ggplot(data_analysis, aes(x = DrinksPerWeek)) + geom_histogram()
p3 <- ggplot(data_analysis, aes(x = SmokeCigs, fill = SmokeCigs)) + geom_bar() + scale_fill_met_d("Derain") 
p4 <- ggplot(data_analysis, aes(x = Age)) + geom_histogram()
p5 <- ggplot(data_analysis, aes(x = HDL_Cholesterol)) + geom_histogram()
p6 <- ggplot(data_analysis, aes(x = Education, fill = Education)) + geom_bar() + scale_fill_met_d("Tam") 
(p1 + p2 + p3)/
  (p4 + p5 + p6)


