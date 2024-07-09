library(gam)

# backfitting algorithm
backfit <- function(x, y, w, s = NULL, tol = 1e-4, maxit = 100) {
  n <- length(y)
  nvar <- ncol(x)
  if(is.null(s)){
    s <- matrix(rep(0, n*nvar), nrow = n, ncol = nvar) #initial smooth components
  }
  r <- y - apply(s, 1, sum)
  fit <- list(fitted.values = 0)
  niter <- 1
  crit <- tol + 1.
  if(is.null(w)) w <- rep(1, n)
  xx <- cbind(rep(1, n), x) # add intercept
  ndig <- -log10(tol) + 1. # number of digits for reporting
  
  while((crit > tol) & (niter < maxit)){
    z <- r + fit$fitted.values
    fit <- lm.wfit(xx, z, w)
    r <- fit$residuals
    deltaf <- 0.
    for(j in 1:nvar){
      old <- s[,j]
      z <- r + s[,j]
      if (j %% 2){
        fit.call <- gam.s(x[,j], z, spar = 1, df = 4)
      }else{
        fit.call <- gam.lo(x[,j], z, span = 0.5, degree = 1, ncols = 1)
      }
      r <- as.double(fit.call$residuals)
      s[,j] <- z - r
      deltaf <- deltaf + weighted.mean((s[,j] - old)^2., w)
    }
    crit <- sqrt(deltaf / sum(w * apply(s, 1., sum)^2.))
    niter <- niter + 1
  }
  
  if(crit > tol) message("Backfitting algorithm failed to converge in ", niter, " iterations.")
  
  fittedvals <- y - r
  return(list("fit" = fit, "fitted.values" = fittedvals, "s" = s, "r" = r))
}

# local scoring algorithm
loc_score <- function(x, y, g, gi, zf, wf, max_iter, tol){
  n <- nrow(x)
  p <- ncol(x)
  alpha <- g(mean(y))
  fold <- rep(0, n)
  fnew <- matrix(0, nrow = n, ncol = p)
  eta <- rep(alpha, n)
  mu <- gi(eta)
  
  delta <- 1
  num_iter <- 0
  while(delta > tol & num_iter < max_iter){
    z <- zf(eta, mu)
    w <- wf(eta, mu)
    
    bf <- backfit(x, z, w, s = NULL)
    fnew <- bf$fitted.values
    eta <- bf$fitted.values
    
    mu <- gi(eta)
    
    delta <- max(abs(fnew - fold))
    fold <- fnew
    num_iter <- num_iter + 1
  }
  
  if(delta > tol) warning("GAM fit failed to converge in ", num_iter, " iterations.")

  return(list(eta = eta, mu = mu, num_iter = num_iter, dev = delta))
}



