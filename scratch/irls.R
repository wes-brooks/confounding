S <- SRE
Lambda <- chol(K)
X <- as.matrix(cbind(1, all.dat[,c('D.Main', 'D.Urb', 'FC', 'MNT', 'MXT', 'Rain')]))

converged <- FALSE
tol <- sqrt(.Machine$double.eps)

dwpr <- glm(resp ~ XX, data=dat, family=poisson(), weights=wt)
eta <- dwpr$linear.predictors
mu <- exp(eta)

M <- rep(0, ncol(Sp))



norm.old <- 0
pseudoCovar <- rbind(cbind(Z, Sp), cbind(Matrix(0, length(M), ncol(X)),  sqrt(tau/2) * diag(length(M))))

while(!converged) {
  z <- eta + (y*exp(-var/2) - mu) / mu
  pseudodata <- c(z, rep(0, length(M)))
  
  lsfit(x=pseudoCovar, y=pseudodata, intercept=FALSE, wt=c(wt*mu * exp(var/2), rep(1, length(M)))) -> f1
  print(f1$coefficients[1:4])
  
  eta <- cbind(X, Sp) %*% f1$coefficients
  mu <- exp(eta)
  
  norm.new <- sum(f1$coefficients^2)
  if (abs(norm.new - norm.old) < tol * (norm.old + tol)) converged <- TRUE
  norm.old <- norm.new
}
