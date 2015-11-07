#------------------------
# Estimate the Cox process model:

cox <- function(y, X, wt, S, Lambda, tau, eta=NULL) {
  converged <- FALSE
  tol <- sqrt(.Machine$double.eps)
  
  if (is.null(eta)) {
    dwpr <- glm(y ~ X, family=poisson(), weights=wt)
    eta <- dwpr$linear.predictors
  }
  mu <- exp(eta)
  u <- rep(0, ncol(S))
  
  norm.old <- 0
  pseudoCovar <- rbind(cbind(X, S %*% Lambda), cbind(Matrix(0, length(u), ncol(X)), sqrt(tau/2) * Diagonal(n=length(u))))
  i <- 0
  
  while(!converged) {
    i <- i+1
    z <- eta + (y - mu) / mu
    pseudodata <- c(z, rep(0, length(u)))
    
    f1 <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, length(u))))
    cat(paste("iter: ", i, '\n', sep=''))
    
    eta <- cbind(X, S %*% Lambda) %*% f1$coefficients
    mu <- exp(eta)
    
    norm.new <- sum(f1$coefficients^2)
    if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
    norm.old <- norm.new
  }
  
  f1$coefficients
}