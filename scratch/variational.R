#-------------------------
# Variational approximation: suppose that u ~ MVN(M, V).
# Here we estimate those parameters:


variational.parameters <- function(y, X, beta, u, S, tau, wt) {
  eta <- as.vector(X %*% beta + S %*% u)
  mu <- exp(eta)
  r <- ncol(S)
  Precision <- tau * diag(nrow=r) + t(S) %*% Diagonal(x=wt * mu) %*% S
  V <- solve(Precision)
  V <- (V + t(V)) / 2
  
  # M <- V %*% (t(S) %*% Diagonal(x=wt) %*% (y - mu) - tau*u) + u
  # M <- E.step(y=y, X=X, beta=beta, S=S, tau=tau, wt=wt)
  M <- rep(0, r)
  
  return(list(M=M, V=V))
}


#------------------------
# Estimate the Cox process model:

E.step <- function(y, X, wt, S, tau, beta) {
  converged <- FALSE
  tol <- sqrt(.Machine$double.eps)
  
  eta <- X %*% beta
  mu <- exp(eta)
  u <- rep(0, ncol(S))
  
  norm.old <- 0
  pseudoCovar <- rbind(S, sqrt(tau/2) * Diagonal(n=length(u)))
  i <- 0
  
  while(!converged) {
    i <- i+1
    z <- eta  - X %*% beta + (y - mu) / mu
    pseudodata <- c(z, rep(0, length(u)))
    
    f1 <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, length(u))))
    cat(paste("iter: ", i, '\n', sep=''))
    
    eta <- X %*% beta + S %*% f1$coefficients
    mu <- exp(eta)
    
    norm.new <- sum(f1$coefficients^2)
    if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
    norm.old <- norm.new
  }
  
  f1$coefficients
}


maximize <- function(y, X, S, wt, M, V) {
  tau <- r/(sum(M^2) + sum(diag(V)))
  
  offset <- as.vector(S %*% M)
  cholV <- t(chol(V))
  pois.exponential <- exp(apply(S, 1, function(z) sum((z %*% cholV)^2) / 2))
  pois.wt <- wt * pois.exponential
  pois.resp <- y / pois.exponential
  gm0 <- glm(pois.resp ~ X - 1, offset=offset, weights=pois.wt, family="poisson")
  
  list(model=gm0, tau=tau)
}

while(diff > .Machine$double.eps) {
  print(i)
  resp <- variational.parameters(y=dat$resp, X=Z, beta=params$model$coefficients, u=rep(0, 209), S=Sp, tau=params$tau, wt=p.wt)
  params <- maximize(y=dat$resp, X=Z, S=Sp, wt=p.wt, M=resp$M, V=resp$V)
  
  diff <- sum((prev - params$model$coefficients)^2)
  print
  prev <- params$model$coefficients
  print(diff)
}

