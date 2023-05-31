#----------------------------------------------------------------------------
# Helper functions
#----------------------------------------------------------------------------


psBasis <- function(x, K = length(x), spline.degree = 3, diff.ord = 2,
                    knots = NULL) {
  if (is.null(knots)) {
    knots.no <- K - spline.degree + 1
    xl <- min(x)
    xr <- max(x)
    xmin <- xl - (xr - xl) / 100
    xmax <- xr + (xr - xl) / 100
    dx <- (xmax - xmin) / (knots.no - 1)
    knots <- seq(xmin - spline.degree * dx, xmax + spline.degree * dx, by = dx)
  }
  X <- splines::spline.des(knots, x, spline.degree + 1, outer.ok = TRUE)$design
  P <- diag(K) # precision
  if (diff.ord > 0) {
    for (d in 1:diff.ord) P <- diff(P)
    P <- crossprod(P)
  }
  return(list(
    X = X, P = P, knots = knots, K = K, spline.degree = spline.degree,
    diff.ord = diff.ord
  ))
}


makeZ <- function(mi){
  z <- matrix(0, nrow = sum(mi), ncol = length(mi))
  k <- 1
  for(i in 1:length(mi)){
    z[k:(k+mi[i]-1),i] <- 1
    k <- k+mi[i]
  }
  z
}


rowrep <- function(X, ntimes){
  #as.matrix(as.data.frame(lapply(as.data.frame(X), rep, ntimes)))
  
  X[rep(seq_along(ntimes), ntimes), ]
}


makeYall <- function(){
  
  X2 <- as.data.frame(rowrep(X, Mi))
  X2$subject <- rep(1:N, Mi)
  X2$subject <- factor(X2$subject)
  Yall <- X2
  Yall$Y <- t(unname(Y))
  class(Yall$Y) = class(Yall$Y)[-1]
  Yall$curve <- 1:nrow(Yall)
  Yall$curve <- factor(Yall$curve)
  Yall
}

fosrcoef <- function(alphaf_post){
  
  S <- ncol(alphaf_post)
    
  alphaf_pci <- apply(alphaf_post[,], 1, function(x) quantile(x, c(.025,.975)))
  alphaf_mean <- rowMeans(alphaf_post[,])

  output <- list(alphaf_pci = alphaf_pci,
                 alphaf_mean = alphaf_mean)
  output
}

fmse <- function(alpha_true, alpha_hat){
  
    sum((alpha_true - alpha_hat)^2)/length(alpha_true)
  
}

mciw <- function(u, l){
  mean(u - l)
  
}

ecp <- function(u, l, alpha){
  
  mean(alpha < u & alpha > l)
  
  
}

fosryhat <- function(m1){
  
  yhat_draws <- array(NA, dim = c(nrow(m1$B), nrow(m1$w_post[[1]]), length(m1$w_post[(S/2 + 1):S])))
  
  for(i in (S/2 + 1):S){
    yhat_draws[,,i - (S/2)] <- m1$B%*%t(rowrep(m1$X,Mi)%*%(m1$alpha_post[[i]][,]) + rowrep(m1$ga_post[[i]], Mi) + m1$w_post[[i]])
  }
  output <- rowMeans( yhat_draws , dims = 2 )
  output
}


