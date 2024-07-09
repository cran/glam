advr <- function(x, y, gi, zf, wf, tol, max_iter){
  n <- nrow(x)
  p <- ncol(x)
  beta_old <- rep(1, p)
  eta_old <- rep(0, n)
  mu_old <- rep(0.5, n)
  dev <- 1
  num_iter <- 0
  while(dev > tol & num_iter < max_iter){
    z <- zf(eta_old, mu_old)
    w <- diag(wf(eta_old, mu_old))
    
    beta_new <- c(solve(t(x) %*% w %*% x) %*% t(x) %*% w %*% z)
    eta_new <- c(x %*% beta_new)
    mu_new <- gi(eta_new)
    mu_new[which(mu_new == 0)] <- 1e-8
    mu_new[which(mu_new == 1)] <- 1 - 1e-8
    
    dev <- sum(abs(beta_new - beta_old))
    
    beta_old <- beta_new
    eta_old <- eta_new
    mu_old <- mu_new
    num_iter <- num_iter + 1
    
  }
  
  if(dev > tol) warning("GLM fit failed to converge in ", num_iter, " iterations.")
  
  return(list(coef = beta_new, eta = eta_new, mu = mu_new, num_iter = num_iter, dev = dev))
}
