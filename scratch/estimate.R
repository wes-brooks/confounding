S.bin <- t(bins) %*% SRE
qrs2 <- qr(S.bin)


loglik <- function(beta, u, ltau, X, y, S, Lambda, wt) {
  eta <- X %*% beta + S %*% Lambda %*% u
  res <- 2 * sum(wt * exp(eta)) - 2 * sum(y * wt * eta) - length(u) * ltau + exp(ltau)*sum(u^2)
  res
}


score.calc <- function(beta, u, ltau, X, y, S, Lambda, wt) {
  eta <- X %*% beta + S %*% Lambda %*% u
  res <- vector()
  
  # Derivative in beta direction:
  res <- c(res, 2*t(X) %*% (wt * (exp(eta) - y)) )
  
  # Derivative in u direction:
  res <- c(res, 2 * t(Lambda) %*% t(S) %*% (wt*(exp(eta) - y)) + 2*exp(ltau)*u)
  
  # Derivative in tau direction:
  # res <- c(res, -r + exp(ltau)*sum(u^2))
  
  res
}


# Approximate log likelihood, marginal to the random effects
marginal.log.lik <- function(beta, u, ltau, X, y, S, Lambda, wt) {
  tau <- exp(ltau)
  eta <- X %*% beta + S %*% Lambda %*% u
  mu <- exp(eta)
  r <- ncol(S)
  half <- Diagonal(x=sqrt(wt * mu)) %*% S %*% Lambda
  to.decompose <- t(half) %*% half + tau * Diagonal(n=r)
  L <- chol(to.decompose)
  ldetL <- determinant(L, logarithm=TRUE)$modulus
  
  
  ll <- sum(y * wt * eta) - sum(wt * mu) - sum(wt*lgamma(y+1)) 
  res <- ll - r/2 * log(2*pi) + r/2*ltau - ldetL - tau/2*sum(u^2)
  res
}


hess <- Matrix(NA, length(par), length(par))
half1 <- Diagonal(x=sqrt(wt * mu)) %*% X
hess[1:2, 1:2] <- t(half1) %*% half1

half2 <- Diagonal(x=sqrt(wt * mu)) %*% S %*% Lambda
hess[3:235, 3:235] <- 2*t(half2) %*% half2 + 2*tau*Diagonal(n=length(u))

hess[1:2, 3:235] <- t(half1) %*% Diagonal(x=sqrt(wt * mu)) %*% S %*% Lambda
hess[3:235, 1:2] <- t(hess[1:2,3:235])


# Setup:
beta <- coef(dwpr)
u <- rep(0, dim(SRE)[1])
Lambda <- chol(K.hat)
S <- SRE.ortho
wt <- p.wt
y <- dat$resp
ltau <- 0
X <- X2

loglik(beta, u, ltau, X, y, S, Lambda, wt)
marginal.log.lik(beta, u, ltau, X, y, S, Lambda, wt)
b <- S %*% Lambda %*% u
Seff <- S %*% Lambda

log.lik <- function(par) {
  beta <- par[1:2]
  u <- par[3:235]
  
  loglik(beta, u, ltau, X, y, S, Lambda, wt)
}


score <- function(par) {
  beta <- par[1:2]
  u <- par[3:235]
  
  score.calc(beta, u, ltau, X, y, S, Lambda, wt)
}


# Standard error formulas: which to use? Square-root whole covariance matrix, or limit to just the fixed effects?
(Diagonal(x=sqrt(wt*mu)) %*% cbind(X, S%*%Lambda) ) %>% qr -> qrX
sqrtm(chol2inv(qrX$qr)[1:2,1:2]) # just fixed effects, than root?
sqrtm(chol2inv(qrX$qr))[1:2,1:2] # or root and then subset?

