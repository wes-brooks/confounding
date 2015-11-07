#--------------------------
# Calculate the (negative) marginal log-likelihod
marginal.laplace <- function(par, y, X, S, u, wt) {
  tryCatch( {
    ltau <- tail(par, 1)
    tau <- exp(ltau)
    p <- ncol(X)
    beta <- par[1:p]
    eta <- X %*% beta + S %*% u
    mu <- exp(eta)
    
    # Calculate the Hessian for the u vector:
    H <- t(S) %*% Diagonal(x=wt * mu) %*% S
    H <- as.matrix((t(H) + H) / 2)
    diag(H) <- diag(H) + tau
    
    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)
    
    r <- ncol(S)
    -sum(wt * (y*eta - exp(eta))) + tau/2 * sum(u^2) - r/2*log(tau) + sum(log(diag(L)))
  }, error=function(e) return(Inf)
  )
}


#-------------------------
# Calculate the gradient of the marginal log-likelihood
gradient.laplace <- function(par, y, X, S, u, wt) {
  tryCatch( {
    ltau <- tail(par, 1)
    tau <- exp(ltau)
    p < -ncol(X)
    beta <- par[1:p]
    eta <- X %*% beta + S %*% u
    mu <- exp(eta)
    r <- ncol(S)
    p <- ncol(X)
    
    # Calculate the Hessian for the u vector:
    H <- t(S) %*% Diagonal(x=wt * mu) %*% S
    H <- as.matrix((t(H) + H) / 2)
    diag(H) <- diag(H) + tau
    
    # Compute the Cholesky decomposition of the Hessian
    L <- chol(H)
    L <- t(L)
    
    deriv <- LogDetDerChol2(L, S, X, mu, wt)
    
    c(as.vector(t(X) %*% Diagonal(x=wt) %*% (y - exp(eta)) + deriv[1:p]), -tau*sum(u^2)/2 + r/2 - tail(deriv, 1))
  }, error=function(e) return(NA)
  )
}




conjugate.gradient.laplace <- function(objective, gradient, y, X, S, beta, u, wt, ltau, verbose=TRUE, tol=sqrt(.Machine$double.eps)) {
  #Initial parameters:
  n = nrow(S)
  r <- ncol(S)
  
  finished <- FALSE
  par <- c(beta, ltau)
  f.new <- objective(par, y, X, S, u, wt)
  f.old <- Inf
  check <- Inf
  
  while(!finished) {
    
    #These iterations restart conjugacy:
    converged = FALSE
    iter <- 0
    while (!converged && iter<100) {
      iter <- iter+1
      
      #Prepare to iterate conjugate gradient descent:
      f.outer <- f.old
      t <- 1
      conv.inner <- FALSE
      i <- 0
      while(f.new < f.old && !converged && !conv.inner && i<p) {
        i <- i+1
        
        dir.new <- gradient(par, y, X, S, u, wt)
        dir.new <- dir.new / sqrt(sum(dir.new^2))
        
        #First iteration, ignore conjugacy - thereafter, use it.
        #par.step is the vector of the new step (in parameter space)
        if (i==1) {  
          par.step <- dir.new
        } else {
          conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
          par.step <- dir.new + conj * s.old
        }
        
        #Find the optimal step size
        #Backtracking: stop when the loss function is majorized
        f.proposed <- objective(par + t*par.step, y, X, S, u, wt)
        condition <- (f.proposed > f.new - sum((t*par.step)*dir.new) - 1/(2*t)*sum((t*par.step)^2))[1]
        while(condition && t > .Machine$double.eps) {
          t = 0.5*t
          f.proposed <- objective(par + t*par.step, y, X, S, u, wt)
          condition = (f.proposed > f.new - sum((t*par.step)*dir.new) - 1/(2*t)*sum((t*par.step)^2))[1]
          
          #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
          if (is.na(condition)) {
            converged = TRUE
            condition = FALSE
          }
        }
        
        #Find the optimal step
        par.proposed <- par + t*par.step
        
        #Make t a little bigger so next iteration has option to make larger step:
        t = t / 0.5 / 0.5
        
        #save for next iteration:
        dir.old <- dir.new
        s.old <- par.step
        
        #Only save the new parameters if they've decreased the loss function
        if (f.proposed < f.old)
          par <- par.proposed
        f.old <- f.new
        f.new <- f.proposed
        
        if ((f.old - f.new) < tol * (tol+abs(f.old))) conv.inner = TRUE
      }
      
      if (verbose) cat(paste("Iteration: ", iter, "; Objective: ", round(f.new, 3), "; Inner iterations: ", i, "\n", sep=""))
      if (abs(f.outer - f.new) < tol * (tol+abs(f.outer))) converged = TRUE
    }
    
    finished <- TRUE
  }
  
  list(par=par, value=f.new)
}




#-------------------------
# Estimate regression parameters for a Cox proess by a Laplace approximation.
cox.laplace <- function(y, X, S, wt, beta.start, tol=sqrt(.Machine$double.eps), verbose=TRUE) {

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
  lik.old <- Inf
  
  # Now iteratively update the estimates of ltau, beta, M, and V until the likelihood converges
  if (verbose) cat("Beginning iteration to estimate beta, u.hat, and tau.\n")
  conv.outer <- FALSE
  while(!conv.outer) {
    
    # Fix ltau and beta, maximize the joint likelihood through u
    if (verbose) cat("Maximizing joint likelihood for u: ")
    converged <- FALSE
    while(!converged) {
      z <- as.vector(eta - X %*% beta + (y - mu) / mu)
      pseudodata <- c(z, rep(0, r))
      pois <- lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu, rep(1, r)))
      
      eta <- as.vector(X %*% beta + S %*% pois$coefficients)
      mu <- exp(eta)
      
      norm.new <- sum(pois$coefficients^2)
      if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
      norm.old <- norm.new
      if (verbose) cat(".")
    }
    u <- pois$coefficients
    if (verbose) cat("done!\n")
    
    # Estimate the regression coefficients
    #res <- optim(beta, fn=marginal.laplace, gr=gradient.laplace, y=dat$resp, X=X, S=S, u=u, ltau=ltau, wt=wt, L=L, method="BFGS")
    #res <- optim(c(beta, ltau), fn=marginal.laplace, gr=gradient.laplace, y=dat$resp, X=X, S=S, u=u, wt=wt, method="BFGS")
    res <- conjugate.gradient.laplace(objective=marginal.laplace, gradient=gradient.laplace, y=dat$resp, X=X, beta=beta, S=S, u=u, ltau=ltau, wt=wt)
    if (verbose) cat(paste("ltau=", round(ltau, 4), "\nbeta=", paste(round(res$par, 3), collapse=', '), '\n', sep=''))
    
    if (verbose) cat(paste("Checking convergence:\n Negative log-likelihood = ", round(res$value, 3), "\n Convergence criterion = ", round(abs(res$value - lik.old) / (tol * (tol + abs(lik.old))), 3), "\n\n"))
    if (abs(res$value - lik.old) < tol * (tol + abs(lik.old))) conv.outer <- TRUE
    beta <- res$par[1:p]
    ltau <- tail(res$par, 1)
    tau <- exp(ltau)
    lik.old <- res$value
    
    eta <- X %*% beta + S %*% u
    mu <- exp(eta)
  }
  
  list(beta=beta, ltau=ltau, u=u)
}