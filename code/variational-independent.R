#-------------------------
# Variatonal lower bound on likelihood (eq 3.1 of Ormerod and Wand, 2012)

likelihood.bound.indep <- function(y, eta, S, wt, ltau, M, diagV) {
  mu <- exp(eta)
  tau <- exp(ltau)
  tryCatch( {
    v <- exp(VariationalVarIndep(diagV, S) / 2)
    
    result <- sum(wt * (y * eta - mu * v))
    result <- result + (r*ltau + sum(log(diagV)) - tau*(sum(M^2) + sum(diagV))) / 2
    return(-result)
  }, error=function(e) return(Inf)
  )
}


newton.indep <- function(y, eta, S, wt, ltau, diagV) {
  tau <- exp(ltau)
  mu <- exp(eta)
  v <- exp(VariationalVarIndep(diagV, S) / 2)
  
  score <- score.indep(mu, S, wt, tau, diagV, v)
  curvature <- curvature.indep(mu, S, wt, diagV, v)
  
  propose <- as.vector(diagV - solve(curvature) %*% score)
  res <- ifelse(propose > 0, propose, diagV)
  
  res
}


#--------------------------
# Gradient of the lower likelihood bound w.r.t. V:
score.indep <- function(mu, S, wt, tau, diagV, v) {
  #grad <- as.vector(colSums(sweep(S^2, 1, wt * mu * v / 2, '*')))
  grad <- -as.vector(t(S^2) %*% as.matrix(wt * mu * v))
  grad <- grad + 1/diagV
  grad <- grad - tau # D_V
  
  grad / 2
}


curvature.indep <- function(mu, S, wt, diagV, v) {
  curve <- -t(S^2) %*% Diagonal(x=wt * mu * v) %*% (S^2) / 4
  diag(curve) <- diag(curve) - diagV^(-2) / 2
  
  curve
}


#----------------------
# Estimate the regression parameters of a Cox process model
# using the variational approximation

cox.variational.indep <- function(y, X, S, wt, beta.start, tol=sqrt(.Machine$double.eps), verbose=TRUE) {
  # Start by estimating an optimal log(tau), assuming the given beta.start and u=rep(0,p)
  beta <- beta.start
  tau <- 100
  r <- ncol(S)
  p <- ncol(X)
  n <- nrow(X)
  
  # Prepare to begin iteration
  pseudoCovar <- rbind(S, sqrt(tau/2) * diag(r))
  eta <- as.vector(X %*% beta)
  mu <- exp(eta)
  norm.old <- Inf
  
  # Iterate to estimate u
  converged <- FALSE
  while(!converged) {
    z <- as.vector(eta - X %*% beta  + (y - mu) / mu)
    pseudodata <- c(z, rep(0, r))
    f1 <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))
    u <- f1$coefficients
    
    eta <- as.vector(X %*% beta + S %*% u)
    mu <- exp(eta)
    
    norm.new <- sum(u^2)
    if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
    norm.old <- norm.new
  }
  
  tau <- r / norm.new
  ltau <- log(tau)
  
  # Initial estimates of M, eta, and mu:
  StS <- apply(S, 1, function(z) sum(z^2))
  M <- u
  eta <- as.vector(X %*% beta + S %*% M)
  mu <- exp(eta)
  
  # Calculate an initial estimate for V
  conv <- FALSE
  ltau.old <- -Inf
  while (!conv) {
    der2 <- -sum(exp(eta + exp(-ltau)*StS/2) * wt* (exp(-ltau)*StS/2 + (exp(-ltau)*StS/2)^2)) - sum(M^2)/2 * exp(ltau)
    der <- sum(exp(eta + exp(-ltau)*StS/2) * wt * exp(-ltau)*StS/2) - exp(ltau)/2 * sum(M^2)
    ltau <- ltau - der/der2
    
    if (abs(ltau - ltau.old) < tol * (tol+abs(ltau.old))) {
      conv <- TRUE
    } else ltau.old <- ltau
  }
  tau <- exp(ltau)
  diagV <- rep(1/tau, r)
  
  # Iterate and maximize the variational approximation to the marginal log-likelihood
  lik.old <- Inf
  conv.outer <- FALSE
  while (!conv.outer) {
    
    # Estimate V (update the variational approximation and maximize it)
    diagV <- newton.indep(y, eta, S, wt, ltau, diagV)
    
    # Estimate M, beta, and tau for fixed V
    v <- exp(VariationalVarIndep(diagV, S) / 2)
    pseudoCovar <- rbind(cbind(X, S), cbind(Matrix(0, length(M), p),  sqrt(tau/2) * diag(r)))
    converged <- FALSE
    while(!converged) {
      z <- eta + (y * v - mu) / mu
      pseudodata <- c(z, rep(0, length(M)))
      pois.model <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v, rep(1, r)))
      
      eta <- cbind(X, S) %*% pois.model$coefficients
      mu <- exp(eta)
      
      norm.new <- sum(pois.model$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new 
    }
    beta <- pois.model$coefficients[1:p]
    M <- tail(pois.model$coefficients, r)
    ltau <- log(r) - log(sum(M^2) + sum(diagV))
    tau <- exp(ltau)
    
    # Check for convergence
    lik <- likelihood.bound.indep(y, eta, S, wt, ltau, M, diagV)
    print(lik.old)
    print(tol)
    print(lik)
    cat(paste("checking conv.outer: (lik - lik.old) / tol = ", round(abs(lik - lik.old) / (tol * (tol + abs(lik.old))), 4), "\n"))
    if (abs(lik - lik.old) < tol * (tol + abs(lik.old))) {
      conv.outer <- TRUE
    } else lik.old <- lik
  }
  
  list(beta=beta, M=M, diagV=diagV, ltau=ltau)
}

