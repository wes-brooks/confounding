# Setting up spatial random effect basis functions
# Use bisquare basis:
bisquare <- function(d, bw)
{
  indx <- which(d < bw)
  resp <- rep(0, length(d))
  resp[indx] <- (1 - (d[indx]/bw)^2)^2

  return(resp)
}


# Set up three different spatial resolutions:
n.1 <- 5
xx <- seq(-0.1, 1.1, len=n.1)
min.1 <- xx[2] - xx[1]
points.1 <- data.frame(x=rep(xx, times=n.1), y=rep(xx, each=n.1))
dist.1 <- sqrt(outer(loc.tot$x, points.1$x, '-')^2 + outer(loc.tot$y, points.1$y, '-')^2)
SRE <- apply(dist.1, 2, function(z) bisquare(z, 1.5*min.1))


n.2 <- round(1.5 * n.1, 0)
xx <- seq(-0.2, 1.2, len=n.2)
min.2 <- xx[2] - xx[1]
points.2 <- data.frame(x=rep(xx, times=n.2), y=rep(xx, each=n.2))
dist.2 <- sqrt(outer(loc.tot$x, points.2$x, '-')^2 + outer(loc.tot$y, points.2$y, '-')^2)
SRE <- cbind(SRE, apply(dist.2, 2, function(z) bisquare(z, 1.5*min.2)))


n.3 <- round(1.5*n.2, 0)
xx <- seq(0.1, 0.9, len=n.3)
min.3 <- xx[2] - xx[1]
points.3 <- data.frame(x=rep(xx, times=n.3), y=rep(xx, each=n.3))
dist.3 <- sqrt(outer(loc.tot$x, points.3$x, '-')^2 + outer(loc.tot$y, points.3$y, '-')^2)
SRE <- cbind(SRE, apply(dist.3, 2, function(z) bisquare(z, 1.5*min.3)))

# Remove any random effect components that are orthgonal to the points
indx <- which(colSums(SRE)==0)
if (length(indx)>0) SRE <- SRE[,-indx]



# Resolution 4 is used for the grid to estimate K and tau
n.4 <- round(2*n.3, 0)
xx <- seq(0, 1, len=n.4)
min.4 <- xx[2] - xx[1]
points.4 <- data.frame(x=rep(xx, times=n.4), y=rep(xx, each=n.4))
dist.4 <- sqrt(outer(loc.tot$x, points.4$x, '-')^2 + outer(loc.tot$y, points.4$y, '-')^2)
bins <- apply(dist.4, 2, function(z) ifelse(z < 1.5*min.4, 1, 0))

# Remove any empty bins
indx <- which(colSums(bins)==0)
if (length(indx)>0) bins <- bins[,-indx]
M <- ncol(bins)

# Compute the method of moments variance Sigma.hat:
D <- Matrix(0, n.tot, M)
for (i in 1:M) D[bins[,i]==1,i] <- residuals(dwpr)[bins[,i]==1]
D2 <- colSums(D) / colSums(bins)
Sigma.hat <- D2 %*% t(D2)
diag(Sigma.hat) <- apply(D, 2, function(z) sum(z^2)) / colSums(bins)


# I _think_ we want the weights here...
#V.bar <- Matrix(diag(sapply(1:M, function(i) sum(1/p.wt * bins[,i]^2) / sum(bins[,i]^2))))
V.bar <- Matrix(diag(M))

S.bar <- t(bins) %*% SRE
qrs <- qr(S.bar)
qrs.R <- qr.R(qrs)
qrs.Q <- qr.Q(qrs)
QQt <- qrs.Q %*% t(qrs.Q)


P <- function(A, QQt) QQt %*% A %*% QQt
Sigma.bar <- P(Sigma.hat, QQt)

# Don't use this for now, since it is resulting in not positive definite K (why?). Just set precision tau to some large number.
# tau <- 1 / lsfit(y=as.vector(Sigma.hat - P(Sigma.hat, QQt)), x=as.vector(V.bar - P(V.bar, QQt)), intercept=FALSE)$coefficients
tau <- 1e9
K.hat <- solve(qrs.R) %*% t(qrs.Q) %*% (Sigma.hat - V.bar/tau) %*% qrs.Q %*% t(solve(qrs.R))
#K.hat <- solve(qrs.R) %*% t(qrs.Q) %*% (Sigma.hat) %*% qrs.Q %*% t(solve(qrs.R))

# Sigma.bar <- SRE %*% K.hat %*% t(SRE) + Matrix(diag(1/p.wt)) / tau
 
# Inverse of the variance-covariance matrix:
V.inv <- Diagonal(x=tau / p.wt)
Precision <- V.inv %*% (Diagonal(n.tot) - SRE %*% solve(solve(K.hat) + t(SRE) %*% V.inv %*% SRE) %*% t(SRE) %*% V.inv)

qrs <- qr(SRE.ortho)
RtR <- t(qr.R(qrs)) %*% qr.R(qrs)
A <- (t(SRE.ortho) %*% V.inv)
Precision <- A %*% SRE - A %*% SRE %*% solve(solve(K.hat) + t(SRE) %*% V.inv %*% SRE) %*% t(SRE) %*% t(A)

# Assume that V = Identity * tau^(-1). Then A %*% SRE is RtR * tau
Precision <- RtR*tau - tau*RtR %*% solve(solve(K.hat) + tau * RtR) %*% RtR * tau


Q <- as.matrix(t(SRE.ortho) %*% Precision %*% SRE.ortho)

