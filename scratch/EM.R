gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

EM <- function(y, X, wt, beta.start, S, Lambda) {
  tol <- sqrt(.Machine$double.eps)
  diff <- Inf
  r <- ncol(S)
  tau <- 1
  u <- rep(0, r)
  beta <- beta.start
  tau.history <- c(tau, tau)
  par.old <- c(beta, tau, u)
  iter <- 0
  while (diff > tol) {
    #if (iter %% 2 ==0) tau <- gm_mean(tail(tau.history, 2))
    iter <- iter + 1
    
    u <- E.step(y, X, wt, S, Lambda, tau, beta)
    beta <- M.step(y, X, wt, S, Lambda, u)
    tau <- length(u) / sum(u^2)
    
    diff <- sum((par.old - c(beta, tau, u))^2)
    par.old <- c(beta, tau, u)
    cat(paste("beta: ", paste(round(beta,3), collapse=', '), "; tau: ", round(tau,3), "\n", sep=''))
    tau.history <- c(tau.history, tau)
  }
  
  par.old
}
# 
# diff <- Inf
# r <- ncol(S)
# tau <- 1
# u <- rep(0, r)
# beta <- dwpr$coefficients
# tau.history <- c(tau, tau)
# par.old <- c(beta, tau, u)
# iter <- 0
# while (diff > tol) {
#   if (iter %% 2 == 0) tau <- gm_mean(tail(tau.history, 2))
#   iter <- iter + 1
#   
#   u <- E.step(dat$resp, Z, p.wt, S, Lambda, tau, beta)
#   beta <- M.step(dat$resp, Z, p.wt, S, Lambda, u)
#   tau <- length(u) / 2 * sum(u^2)
#   
#   diff <- sum((par.old - c(beta, tau, u))^2)
#   par.old <- c(beta, tau, u)
#   cat(paste("beta: ", paste(round(beta,3), collapse=', '), "; tau: ", round(tau,3), "\n", sep=''))
#   tau.history <- c(tau.history, tau)
# }
# par.non.orthogonal <- par.old
