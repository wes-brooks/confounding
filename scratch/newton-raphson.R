#-------------------------
# Variational approximation: suppose that u ~ MVN(M, V).
# Estimate parameters by Newton-Raphson iteration:

library(matrixcalc)

#--------------------------------
# Utility functions:
Q <- function(A) {
  p <- ncol(A)
  
  #Element-wise multiplication of kronecker product matrices:
  kronecker(A, matrix(1, 1, p)) * kronecker(matrix(1, 1, p), A)
}


#--------------------------
# Calculate the Hessian of the parameters:
hessian <- function(X, S, b, M, V) {
  p <- ncol(X)
  r <- ncol(S)
  
  eta <- as.vector(X %*% b + S %*% M)
  half <- t(chol(V))
  B <- Diagonal(x=eta + 0.5 * apply(S, 1, function(z) sum((z %*% half)^2)))
  Q.S <- Q(S)
  D.r <- D.matrix(r)
  
  # Calculate the Hessian in blocks:
  H.b.b <- t(X) %*% B %*% X
  
  H.b.M <- -t(X) %*% B %*% S
  H.b.V <- -0.5 * t(X) %*% B %*% Q.S 
  
}
 
