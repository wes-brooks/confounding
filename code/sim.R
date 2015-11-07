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
n.q <- 20000
inc <- rexp(n.q + 1)
x <- cumsum(inc[1:n.q]) / sum(inc)
y <- runif(n.q)
loc <- data.frame(x=x, y=y)

# Generate the quadrature points
n <- round(sqrt(10000))
x <- rep(seq(0, 1, length=n), each=n)
y <- rep(seq(0, 1, length=n), times=n)
loc <- rbind(loc, data.frame(x=x, y=y))


#--------------------
# Simulate the random effect
re <- GRF(loc=data.frame(x=loc$x, y=loc$y), tau=0.1, sigma=0.5)




#---------------------
# Simulate the covariates
p = 3
X = matrix(0, nrow(loc), 0)
fields <- list()
for(j in 1:p) {
  d <- GRF(loc=data.frame(x=loc$x, y=loc$y), tau=0.25, sigma=1)
  fields[[j]] <- d
  X = cbind(X, d@data[[1]])
}

# X = as.matrix((x + y) / 4)

# 
# 
# #--------------------
# # Apply correlation between covariates, random effect:
# rho <- matrix(0.4, 2, 2)
# diag(rho) <- 1
# Z <- cbind(X, re@data[[1]]) %*% chol(rho)
# 




obs.indx <- 1:n.q
quad.indx <- (n.q+1):nrow(loc)
# eta <- 2 + Z[obs.indx,1] + Z[obs.indx,2] #+ X[,1]
eta <- 5 + 0.2*X[,1] + 0.1*X[,2]  + re@data[[1]]
sup <- exp(max(eta))


#-----------------------
# Apply the simulated regression model and then resample the points to get the Cox process:

count <- 0
while(TRUE) {
  count <- count + 1
  

  
  # Grab some points for the homogeneous Poisson process
  n.ppp <- rpois(1, sup)
  indx.ppp <- sample(obs.indx, n.ppp)
  X.ppp <- X[indx.ppp, ]
  
  write.csv(loc[indx.ppp,], file=paste("~/Dropbox/confounding/output/loc.ppp", count, "out", sep="."), row.names=FALSE)
  
  # Thin the points we just selected
  prob <- exp(eta[indx.ppp]) / sup
  cox.indx <- which(rbinom(n.ppp, size=1, prob=prob)==1)
  X.cox <- as.matrix(X.ppp[cox.indx,])
  n.cox <- length(cox.indx)
  loc.cox <- loc[indx.ppp[cox.indx],]
  
  write.csv(loc[indx.ppp[cox.indx],], file=paste("~/Dropbox/confounding/output/loc.cox", count, "out", sep="."), row.names=FALSE)
  
  # Get the locations remaining after the homogeneous Poisson process
  X.q <- X[quad.indx,]
  n.q <- length(quad.indx)
  loc.q <- loc[quad.indx,]
  
  n.tot <- n.q + n.cox
  loc.tot <- rbind(loc.cox, loc.q)
  
  # Estimate the model by DWPR
  presence <- c(rep(1, n.cox), rep(0, n.q))
  #dd.tot <- deldir(loc.tot)
  #p.wt <- dd.tot$summary$dir.area
  p.wt <- rep(1/n.tot, n.tot)
  XX <- rbind(X.cox, X.q)
  dat <- data.frame(resp=round(presence/p.wt, 0), XX=XX, x=loc.tot$x, y=loc.tot$y, wt=p.wt)
  rownames(dat) <- 1:n.tot
  dat$id <- rownames(dat)
  
  dwpr <- glm(resp~XX, family=poisson(), weights=wt, data=dat)
  #  summary(dwpr)
  
  sink(file=paste("~/Dropbox/confounding/output/dwpr", count, "out", sep="."), append=FALSE)
  print(dwpr$coefficients)
  sink()


  Z <- as.matrix(cbind(1, dat[,c('XX.1', 'XX.2', 'XX.3')]))
  
  
  srecov <- SRE.covariance(sre=spatial, dwpr=dwpr, loc=loc.tot)
  Lambda <- srecov$Lambda
  S <- srecov$S
  
  #-----------------------
  # Now estimate the model with random effects orthogonal to the fixed effects:
  srecov.o <- SRE.covariance(sre=spatial, dwpr=dwpr, loc=loc.tot, orthogonal=TRUE, X=as.data.frame(Z))
  Lambda.o <- srecov.o$Lambda
  S.o <- srecov.o$S
  

  ortho <- EM(y=dat$resp, X=Z, wt=p.wt, S=S.o, Lambda=Lambda.o, beta.start=dwpr$coefficients)
  ortho <- EM(y=dat$resp, X=Z, wt=p.wt, Spatial=SL2, beta.start=dwpr$coefficients)
  
  sink(file=paste("~/Dropbox/confounding/output/ortho", count, "out", sep="."), append=FALSE)
  cat(ortho)
  sink()
  
#   marg.o <- function(ltau) {
#     m1.o <- cox(y=dat$resp, X=Z, wt=p.wt, loc=loc.tot, Lambda=Lambda.o, S=S.o, tau=exp(ltau))
#     beta.o <- m1.o[1:4]
#     u.o <- m1.o[5:length(m1.o)]
#     
#   
#   
#     res <- marginal.log.lik(beta.o, u.o, ltau, Z, y, S.o, Lambda.o, p.wt)
#     cat("ltau: ", ltau, "; ll: ", res, '\n', sep='')
#     res
#   }
  
  
  
  #-----------------------
  # Now estimate the model with random effects not orthogonal to the fixed effects:
  
  #-----------------------
  # Estimate the matrix K
  qrs <- qr(t(bins) %*% S / colSums(bins))
  qrs.R.inv <- solve(qr.R(qrs))
  qrs.Q <- qr.Q(qrs)
  K <- qrs.R.inv %*% t(qrs.Q) %*% Sigma.hat %*% qrs.Q %*% t(qrs.R.inv)
  Lambda <- chol(K)
  
  
  non.ortho <- EM(y=dat$resp, X=Z, wt=p.wt, S=S, Lambda=Lambda, beta.start=dwpr$coefficients)
  
  
  sink(file=paste("~/Dropbox/confounding/output/non.ortho", count, "out", sep="."), append=FALSE)
  cat(non.ortho)
  sink()
  
  # eta <- dwpr$linear.predictors
  
#   marg <- function(ltau) {
#     m1 <- cox(y=dat$resp, X=Z, wt=p.wt, loc=loc.tot, Lambda=Lambda, S=S, tau=exp(ltau))
#     beta <- m1[1:4]
#     u <- m1[5:length(m1)]
#     
#     
#     
#     res <- marginal.log.lik(beta, u, ltau, Z, y, S, Lambda, p.wt)
#     cat("ltau: ", ltau, "; ll: ", res, '\n', sep='')
#     res
#   }
# 
}
