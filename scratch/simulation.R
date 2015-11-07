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


# #------------------
# # Parameters
# intercept <- c(1,2,5,10)
# beta <- c(1, -1)
# rando.corr <- c(0, 0)
# rho <- matrix(0, p, p)
# rho[1,] <- rho[,1] <- 0.2
# diag(rho) <- 1


#--------------------
# Generate the observation points
n <- 10000
inc <- rexp(n + 1)
x <- cumsum(inc[1:n]) / sum(inc)
y <- runif(n)
loc <- data.frame(x=x, y=y)

# Generate the quadrature points
n <- round(sqrt(10000))
inc <- rexp(n + 1)
x <- rep(seq(0, 1, length=n), each=n)
y <- rep(seq(0, 1, length=n), times=n)
loc <- rbind(loc, data.frame(x=x, y=y))


#--------------------
# Simulate the random effect
re <- GRF(loc=data.frame(x=loc$x, y=loc$y), tau=0.15, sigma=1)




#---------------------
# Simulate the covariates
p = 1
X = matrix(0, nrow(loc), 0)
fields <- list()
for(j in 1:p) {
  d <- GRF(loc=data.frame(x=loc$x, y=loc$y))
  fields[[j]] <- d
  X = cbind(X, d@data[[1]])
}

# X = as.matrix((x + y) / 4)

# 
# 
# #--------------------
# # Apply correlation between covariates, random effect:
rho <- matrix(0.4, 2, 2)
diag(rho) <- 1
Z <- cbind(X, re@data[[1]]) %*% chol(rho)





#-----------------------
# Apply the simulated regression model and then resample the points to get the Cox process:
# for (j in 2:p)
# {
obs.indx <- 1:10000
  eta <- 2 + Z[obs.indx,1] + Z[obs.indx,2] #+ X[,1]
  sup <- exp(max(eta))
  
  # Grab some points for the homogeneous Poisson process
  n.ppp <- rpois(1, sup)
  indx.ppp <- sample(obs.indx, n.ppp)
  X.ppp <- Z[indx.ppp, j]
  
  # Thin the points we just selected
  prob <- exp(eta[indx.ppp]) / sup
  cox.indx <- which(rbinom(n.ppp, size=1, prob=prob)==1)
  X.cox <- as.matrix(X.ppp[cox.indx])
  n.cox <- length(cox.indx)
  loc.cox <- loc[indx.ppp[cox.indx],]
  
  # Get the locations remaining after the homogeneous Poisson process
  X.q <- Z[-obs.indx,j]
  n.q <- length(X.q)
  loc.q <- loc[-obs.indx,]
  
  n.tot <- n.q + n.cox
  loc.tot <- rbind(loc.cox, loc.q)
  
  # Estimate the model by DWPR
  presence <- c(rep(1, n.cox), rep(0, n.q))
  #dd.tot <- deldir(loc.tot)
  #p.wt <- dd.tot$summary$dir.area
  p.wt <- rep(1/n.tot, n.tot)
  XX <- c(X.cox, X.q)
  dat <- data.frame(resp=round(presence/p.wt, 0), XX=XX, x=loc.tot$x, y=loc.tot$y, wt=p.wt)
  rownames(dat) <- 1:n.tot
  dat$id <- rownames(dat)
  
  dwpr <- glm(resp~XX, family=poisson(), weights=wt, data=dat)
  summary(dwpr)
#  }


#----------------------
# Put the species points and quadrature points in a Voronoi diagram
D.tot <- Matrix(0, n.tot, n.tot)
for (i in 1:(n.tot-1)) {indx <- dd.tot$delsgs$ind1[dd.tot$delsgs$ind2==i]; if (length(indx) > 0) D.tot[i,indx] <- 1}
off.diag <- D.tot + t(D.tot)
diag.part = Matrix(0, n.tot, n.tot)
diag(diag.part) <- off.diag %*% matrix(1, n.tot, 1)
Q <- diag.part - off.diag

begin <- Sys.time()
chol.D <- chol(D)
Sys.time() - begin

begin <- Sys.time()
X <- cbind(1, X)
P.c <- Matrix(diag(n)) - X %*% solve(t(X) %*% X) %*% t(X)
Sys.time() - begin

diag(D) <- D %*% matrix(1, n.tot)




