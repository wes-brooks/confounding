# Replicating some analyses from the Renner et al. paper
#-------------------

# Import libraries:
library(dplyr)
library(spatstat)
library(ppmlasso)
library(RandomFields)
library(spatstat)
library(INLA)



# Set up the study window:
load("data/ppm/Eucalyptus sparsifolia Atlas 2012.RData") #Contains X and Y
load("data/ppm/Quad1000.RData") #Contains quad
ux = sort(unique(quad$X))
uy = sort(unique(quad$Y))
nx = length(ux)
ny = length(uy)
col.ref = match(quad$X, ux)
row.ref = match(quad$Y, uy)
all.vec = rep(NA, max(row.ref)*max(col.ref))
vec.ref = (col.ref - 1)*max(row.ref) + row.ref
all.vec[vec.ref] = 1
Sydney.mask = matrix(all.vec, max(row.ref), max(col.ref), dimnames=list(uy, ux))
Sydney.win = as.owin(im(Sydney.mask, xcol=ux, yrow=uy))


# Make point pattern and quadrature scheme:
ppp.dat = ppp(X, Y, window = Sydney.win, check = FALSE)
quads = ppp(quad$X, quad$Y, window = Sydney.win)
Q = quadscheme(data=ppp.dat, dummy=quads, method="grid", ntile=c(nx, ny), npix=c(nx, ny))


# Set up covariate lists:
X.des = cbind(poly(quad$FC, quad$MNT, quad$MXT, quad$Rain, degree = 2,
                   raw = TRUE), poly(sqrt(quad$D.Main), sqrt(quad$D.Urb), degree = 2,
                                     raw = TRUE), quad$soil)
int.list = list()
for (i in 1:dim(X.des)[2]) {
  all.vec = rep(NA, max(row.ref)*max(col.ref))
  vec.ref = (col.ref - 1)*max(row.ref) + row.ref
  all.vec[vec.ref] = X.des[,i]
  int.list[[i]] = im(matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux)), xcol = ux, yrow = uy)
}

names(int.list) = paste("V", 1:dim(X.des)[2], sep = "")
pred.list = int.list
set.0 = 15:19 #Variables to set to 0
for (v in set.0) {
  pred.list[[v]]$v = 0*pred.list[[v]]$v
}



# Establish observation and quadrature points:
sp.xy <- cbind(X,Y)
quad.xy <- as.matrix(quad[,c('X','Y')])
data.xy <- rbind(sp.xy, quad.xy)


# Matrices of data:
sp.dat = data.frame(X, Y, D.Main, D.Urb, FC, MNT, MXT, Rain)
quad.dat = quad[,c('X', 'Y', 'D.Main', 'D.Urb', 'FC', 'MNT', 'MXT', 'Rain')]
sp.dat$Pres = 1
quad.dat$Pres = 0
all.dat = data.frame(rbind(sp.dat, quad.dat))


# Now estimate the model parameters:
p.wt <- rep(1/nrow(all.dat), nrow(all.dat))
p.wt[all.dat$Pres==0] <- 86227 / sum(all.dat$Pres==0)
all.dat$resp <- round(all.dat$Pres/p.wt)
dwpr <- glm(resp ~ D.Main + D.Urb + FC + MNT + MXT + Rain, weights=p.wt, data=all.dat, family=poisson())
summary(dwpr)






# Constructing the mesh and defining the SPDE
hull = inla.nonconvex.hull(as.matrix(data.xy), convex = -0.01, resolution = 500)
hull <- Polygon(hull$loc)
colnames(hull@coords) <- c('x', 'y')
XX <- range(all.dat$X)
YY <- range(all.dat$Y)


#-------------------------
# Set up a spatial random effect model:
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
n.1 <- 25
len.xx <- (n.1 / (diff(YY) / diff(XX))) %>% sqrt %>% round
len.yy <- (n.1 / (diff(XX) / diff(YY))) %>% sqrt %>% round
xx <- seq(min(all.dat$X), max(all.dat$X), len=len.xx)
yy <- seq(min(all.dat$Y), max(all.dat$Y), len=len.yy)
min.1 <- min(xx[2] - xx[1], yy[2] - yy[1])
points.1 <- data.frame(x=rep(xx, times=len.yy), y=rep(yy, each=len.xx))
dist.1 <- sqrt(outer(all.dat$X, points.1$x, '-')^2 + outer(all.dat$Y, points.1$y, '-')^2)
SRE <- apply(dist.1, 2, function(z) bisquare(z, 1.5*min.1))


n.2 <- round(3 * n.1)
len.xx <- (n.2 / (diff(YY) / diff(XX))) %>% sqrt %>% round
len.yy <- (n.2 / (diff(XX) / diff(YY))) %>% sqrt %>% round
xlim <- c(0.05, 0.95) * diff(XX) + min(XX)
ylim <- c(0.05, 0.95) * diff(YY) + min(YY)
xx <- seq(xlim[1], xlim[2], len=len.xx)
yy <- seq(ylim[1], ylim[2], len=len.yy)
min.2 <- min(xx[2] - xx[1], yy[2] - yy[1])
points.2 <- data.frame(x=rep(xx, times=len.yy), y=rep(yy, each=len.xx))
dist.2 <- sqrt(outer(all.dat$X, points.2$x, '-')^2 + outer(all.dat$Y, points.2$y, '-')^2)
SRE <- cbind(SRE, apply(dist.2, 2, function(z) bisquare(z, 1.5*min.2)))


n.3 <- round(3 * n.2)
len.xx <- (n.3 / (diff(YY) / diff(XX))) %>% sqrt %>% round
len.yy <- (n.3 / (diff(XX) / diff(YY))) %>% sqrt %>% round
xlim <- c(0.02, 0.98) * diff(XX) + min(XX)
ylim <- c(0.02, 0.98) * diff(YY) + min(YY)
xx <- seq(xlim[1], xlim[2], len=len.xx)
yy <- seq(ylim[1], ylim[2], len=len.yy)
min.3 <- min(xx[2] - xx[1], yy[2] - yy[1])
points.3 <- data.frame(x=rep(xx, times=len.yy), y=rep(yy, each=len.xx))
dist.3 <- sqrt(outer(all.dat$X, points.3$x, '-')^2 + outer(all.dat$Y, points.3$y, '-')^2)
SRE <- cbind(SRE, apply(dist.3, 2, function(z) bisquare(z, 1.5*min.3)))

# Remove any random effect components that are orthgonal to the points
indx <- which(colSums(SRE)==0)
if (length(indx)>0) SRE <- SRE[,-indx]



# Resolution 4 is used for the grid to estimate K and tau
n.4 <- round(4 * n.3)
len.xx <- (n.4 / (diff(YY) / diff(XX))) %>% sqrt %>% round
len.yy <- (n.4 / (diff(XX) / diff(YY))) %>% sqrt %>% round
xlim <- XX
ylim <- YY
xx <- seq(xlim[1], xlim[2], len=len.xx)
yy <- seq(ylim[1], ylim[2], len=len.yy)
min.4 <- min(xx[2] - xx[1], yy[2] - yy[1])
points.4 <- data.frame(x=rep(xx, times=len.yy), y=rep(yy, each=len.xx))
dist.4 <- sqrt(outer(all.dat$X, points.4$x, '-')^2 + outer(all.dat$Y, points.4$y, '-')^2)
bins <- apply(dist.4, 2, function(z) ifelse(z < 1.5*min.4, 1, 0))

# Remove any empty bins
indx <- which(colSums(bins)==0)
if (length(indx)>0) bins <- bins[,-indx]
M <- ncol(bins)


# Compute the method of moments variance Sigma.hat:
D <- Matrix(0, nrow(all.dat), M)
for (i in 1:M) D[bins[,i]==1,i] <- residuals(dwpr)[bins[,i]==1]
D2 <- colSums(D) / colSums(bins)
Sigma.hat <- D2 %*% t(D2)
diag(Sigma.hat) <- apply(D, 2, function(z) sum(z^2)) / colSums(bins)

# Estimate the matrix K
S.bar <- t(bins) %*% SRE / colSums(bins)
qrs <- qr(S.bar)
qrs.R <- qr.R(qrs)
qrs.Q <- qr.Q(qrs)
QQt <- qrs.Q %*% t(qrs.Q)


P <- function(A, QQt) QQt %*% A %*% QQt

K <- solve(qrs.R) %*% t(qrs.Q) %*% Sigma.hat %*% qrs.Q %*% t(solve(qrs.R))



# Standard error formulas: which to use? Square-root whole covariance matrix, or limit to just the fixed effects?
(Diagonal(x=sqrt(p.wt*mu)) %*% cbind(X, S%*%Lambda) ) %>% qr -> qrX
sqrtm(chol2inv(qrX$qr)[1:ncol(X),1:ncol(X)]) # just fixed effects, than root?
