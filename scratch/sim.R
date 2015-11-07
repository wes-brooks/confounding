# Simulate a Gaussian random field a set of points

# Import libraries
library(RandomFields)
library(expm)
library(ggplot2)
library(deldir)
library(dplyr)

#---------------
# Function to simulate GRFs:
GRF <- function(loc, sigma=1, tau=0.1, nugget=0, fixed=null) {
  covariate.model <- RMexp(var=sigma^2, scale=tau)
  if (nugget!=0) 
    covariate.model <- covariate.model + RMnugget(var=nugget)
  
  dat <- RFsimulate(covariate.model, x=loc$x, y=loc$y)
  dat
}
 

#---------------------
# Generate points from a homogeneous Poisson process:
lam <- 4000
A <- (1-0) * (1 -0)
n <- rpois(1, A * lam)
inc.x <- rexp(n+1)
x <- cumsum(inc.x[1:n]) / sum(inc.x)
y <- runif(n)
loc <- data.frame(x=x, y=y)


#---------------------
# Simulate the covariates
p = 1
X = matrix(0, n, 0)
fields <- list()
for(j in 1:p) {
  d <- GRF(loc=data.frame(x=x, y=y))
  fields[[j]] <- d
  X = cbind(X, d@data[[1]])
}



#--------------------
# Simulate the random effect
re <- GRF(loc=data.frame(x=x, y=y))



#--------------------
# Apply correlation between covariates, random effect:
rando.corr <- c(0, 0.8, 0)
rho <- matrix(0.4, p, p)
rho <- rbind(rho, rando.corr)
rho <- cbind(rho, c(rando.corr, 1))
diag(rho) <- 1
Z <- cbind(X, re@data[[1]]) %*% chol(rho)
rando <- Z[,4]
Z <- Z[,1:3]


#-----------------------
# Apply the simulated regression model and then resample the points to get the Cox process:
alpha <- 1
#beta <- as.matrix(c(0.5, -0.75, -0.2))
beta <- as.matrix(2)
#eta <- alpha + Z %*% beta + rando
eta <- alpha + X %*% beta
lambda <- exp(eta) / exp(max(eta))
cox.indx <- which(rbinom(n, size=1, prob=lambda)==1)
#Z.cox <- Z[cox.indx,]
X.cox <- as.matrix(X[cox.indx,])
n.cox = length(cox.indx)
loc.cox <- loc[cox.indx,]


#--------------------
# Generate the quadrature points
n.q <- 100000
inc.q <- rexp(n.q+1)
x.q <- cumsum(inc.q[1:n.q]) / sum(inc.q)
y.q <- runif(n.q)
loc.q <- data.frame(x=x.q, y=y.q)


#----------------------
# Put the species points and quadrature points in a Voronoi diagram
loc.tot <- rbind(loc.cox, loc.q)
n.tot <- nrow(loc.tot)
dd.tot <- deldir(loc.tot)
D.tot <- Matrix(0, n.tot, n.tot)
for (i in 1:n.tot) {indx <- dd.tot$delsgs$ind1[dd.tot$delsgs$ind2==i]; if (length(indx) > 0) D.tot[i,indx] <- 1}
D <- D.tot + t(D.tot)
diag(D) <- D %*% matrix(1, n.tot)


#---------------------
# Simulate the covariates
X.q = matrix(0, n.q, 0)
for(j in 1:p) {
  d <- RFsimulate(covariate.model, x=x.q, y=y.q, data=fields[[j]])
  X.q = cbind(X.q, d@data[[1]])
}
Z.q <- X.q %*% chol(rho[1:p,1:p])


#-------------
# Estimate the model by DWPR
presence <- c(rep(1, n.cox), rep(0, n.q))
#p.wt <- dd.tot$summary$dir.area
p.wt <- c(rep(1e-6, n.cox), rep(1/n.tot, n.q))
#ZZ = rbind(Z.cox, Z.q)
#dwpr <- glm(presence/p.wt ~ ZZ, family=poisson(), weights=p.wt)
XX <- rbind(X.cox, X.q)
XX <- sweep(XX, 2, apply(XX, 2, mean)) %>% (function(z) sweep(z, 2, apply(z, 2, sd), '/'))
dwpr <- glm(round(presence/p.wt,0) ~ XX, family=poisson(), weights=p.wt)
summary(dwpr)
