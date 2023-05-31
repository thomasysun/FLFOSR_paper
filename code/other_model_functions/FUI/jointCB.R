## Compare the computing time of two methods to obtain joint confidence bands

rm(list = ls())
library(mvtnorm)
library(ggplot2)

B <- 300 ## number of samples to simulate

m.list <- seq(100, 2500, 100)
comp.time <- matrix(NA, nrow = length(m.list), ncol = 3)
colnames(comp.time) <- c("m", "traditional", "fpca")
comp.time[,1] <- m.list
for(i in 1:length(m.list)){
  m <- m.list[i]
  
  ## simulate samples 
  est <- matrix(rnorm(B*m), nrow = B, ncol = m)
  est <- fpca.face(est)$Yhat ## smooth the simulated estimates
  Sigma <- var(est) ## variance of the simulated data
  
  N <- 10000 ## sample size in simulation-based approach
  ## method 1: Traditional
  ptm1 <- proc.time()
  x.sample1 <- rmvnorm(N, mean = colMeans(est), sigma = Sigma)
  comp.time[i,2] <- (proc.time() - ptm1)[3]
  
  ## method 2: FPCA-based
  ptm2 <- proc.time()
  fit_fpca <- fpca.face(est)
  ## extract estimated eigenfunctions/eigenvalues
  phi <- fit_fpca$efunctions
  lambda <- fit_fpca$evalues
  K <- length(fit_fpca$evalues)
  ## simulate random coefficients 
  theta <- matrix(rnorm(N*K), nrow=N, ncol=K) # generate independent standard normals 
  theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
  X_new <- theta %*% t(phi) # simulate new functions
  x.sample2 <- X_new + t(fit_fpca$mu %o% rep(1,N)) # add back in the mean function
  comp.time[i,3] <- (proc.time() - ptm2)[3]
  
  print(paste0("m = ", m, ", ", i, "/", length(m.list)))
  save(comp.time, file = "./jointCB.rda")
}


## plot the computing time of two approaches
comp.time <- as.data.frame(comp.time)
comp.time.plot <- comp.time[1:25,]
plt <- ggplot(comp.time.plot) +
  theme_bw() +
  geom_line(aes(x = m, y = traditional, color = "Traditional method"), lwd = 1.2) +
  geom_line(aes(x = m, y = fpca, color = "New method"), lwd = 1.2) +
  scale_color_manual(values = c("blue", "red"), name = "") +
  labs(y = "Time (s)", x = "Dimension of the functional domain (L)", 
       title = "Computing Time for Joint Confidence Bands") +
  theme(legend.position = c(0.13, 0.89),
        legend.background = element_rect(fill=alpha('white', 0)),
        plot.title = element_text(face = "bold", hjust = 0.5))
plt

