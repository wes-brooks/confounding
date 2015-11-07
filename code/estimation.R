#----------------------
# Estimate the regression parameters of a Cox process model
# using the variational approximation

cox.variational <- function(y, X, S, wt, beta.start, tol=sqrt(.Machine$double.eps), verbose=TRUE) {
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
    z <- as.vector(eta - X %*% beta + (y - mu) / mu)
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
  V <- 1/tau * Diagonal(r)
  
  # Iterate and maximize the variational approximation to the marginal log-likelihood
  lik.old <- Inf
  conv.outer <- FALSE
  while (!conv.outer) {
    
    # Estimate V (update the variational approximation and maximize it)
    result <- conjugate.gradient(y, eta, S, wt, ltau, M, V, objective=likelihood.bound, gradient=score, tol=tol)
    V <- result$V
    
    # Estimate M, beta, and tau for fixed V
    if (verbose) cat("Holding V fixed to estimate M, tau, and beta")
    cholV <- as.matrix(t(chol(V)))
    v <- exp(VariationalVar(cholV, S) / 2)
    pseudoCovar <- rbind(cbind(X, S), cbind(Matrix(0, length(M), p),  sqrt(tau/2) * diag(r)))
    converged <- FALSE
    while(!converged) {
      z <- eta + (y * v - mu) / mu
      pseudodata <- c(z, rep(0, length(M)))
      pois.model <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt * mu * v, rep(1, r)))
      
      eta <- cbind(X, S) %*% pois.model$coefficients
      mu <- exp(eta)
      
      if (verbose) cat(".")
      norm.new <- sum(pois.model$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new 
    }
    beta <- pois.model$coefficients[1:p]
    M <- tail(pois.model$coefficients, r)
    ltau <- log(r) - log(sum(M^2) + sum(diag(V)))
    tau <- exp(ltau)
    if (verbose) cat(paste("done!\n beta = ", beta, "\n ltau = ", ltau, "\n", sep=''))
    
    # Check for convergence
    lik <- likelihood.bound(y, eta, S, wt, ltau, M, V)
    cat(paste("Checking final convergence criterion = ", round(abs(lik - lik.old) / (tol * (tol + abs(lik.old))), 4), "\n"))
    if (abs(lik - lik.old) < tol * (tol + abs(lik.old))) {
      conv.outer <- TRUE
    } else lik.old <- lik
  }
  
  list(beta=beta, M=M, V=V, ltau=ltau)
}