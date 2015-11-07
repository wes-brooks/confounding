#------------------------
# Estimate the Cox process model:

E.step <- function(y, X, wt, S, Lambda, tau, beta) {
  converged <- FALSE
  tol <- sqrt(.Machine$double.eps)
  
  eta <- X %*% beta
  mu <- exp(eta)
  u <- rep(0, ncol(S))
  
  norm.old <- 0
  pseudoCovar <- rbind(S %*% Lambda, sqrt(tau/2) * Diagonal(n=length(u)))
  i <- 0
  
  while(!converged) {
    i <- i+1
    z <- eta  - X %*% beta + (y - mu) / mu
    pseudodata <- c(z, rep(0, length(u)))
    
    f1 <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, length(u))))
    cat(paste("iter: ", i, '\n', sep=''))
    
    eta <- X %*% beta + S %*% Lambda %*% f1$coefficients
    mu <- exp(eta)
    
    norm.new <- sum(f1$coefficients^2)
    if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
    norm.old <- norm.new
  }
  
  f1$coefficients
}

M.step <- function(y, X, wt, S, Lambda, u) {
  as.vector(S %*% Lambda %*% u) -> J
  glm(y~X-1, offset=J, weights=wt, family='poisson')$coef-> beta
  beta
}


