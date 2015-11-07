#-------------------------
# Variatonal lower bound on likelihood (eq 3.1 of Ormerod and Wand, 2012)
likelihood.bound <- function(y, eta, S, wt, ltau, M, V) {
  mu <- exp(eta)
  tau <- exp(ltau)
  tryCatch( {
    cholV <- chol(V)
    cholV <- t(as.matrix(cholV))
    v <- exp(VariationalVar(cholV, S) / 2)
    
    result <- sum(wt * (y * eta - mu * v))
    result <- result + (r*ltau + 2*sum(log(diag(cholV))) - tau*(sum(M^2) + sum(diag(V)))) / 2
    return(-result)
  }, error=function(e) return(Inf)
  )
}



#--------------------------
# Gradient of the lower likelihood bound:
score <- function(y, eta, S, wt, ltau, V) {
  r <- ncol(S)
  tau <- exp(ltau)
  mu <- exp(eta)
  
  cholV <- chol(V)
  cholV <- t(as.matrix(cholV))
  v <- exp(VariationalVar(cholV, S) / 2)
  
  grad <- VariationalScore(mu, wt, tau, v, as.matrix(V), S)
  
  as.vector(grad)
}