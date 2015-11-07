# Gamma hyperparameters for tau
a <- 0
b <- 0

# Normal hyperparameters for beta:
mu <- 0
tau <- 0

#Likelihood:
eta <- X2 %*% beta + SRE %*% u

grad.u <- t(SRE) %*% Diagonal(x=p.wt) %*% (dat$resp - exp(eta)) + tau * (Precision %*% u)
grad.beta <- t(X2) %*% Diagonal(x=p.wt) %*% (dat$resp - exp(eta))
grad.tau <- r/2/tau + u %*% Precision %*% u

ll <- function(beta, tau) {
  sum(p.wt*(dat$resp*eta - exp(eta))) - tau/2 * t(u) %*% Precision %*% u + 1/2 * sum(log(eigen(Precision)$values)) + r/2 * log(tau)
}

