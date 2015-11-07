# Replicating some analyses from the Renner et al. paper
#-------------------

# Import libraries:
library(dplyr)
library(spatstat)
library(ppmlasso)
library(RandomFields)
library(spatstat)



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



# Constructing the mesh and defining the SPDE
hull = inla.nonconvex.hull(as.matrix(data.xy), convex = -0.01, resolution = 800)
mesh = inla.mesh.2d(boundary = hull, cutoff = 3, offset = c(20, 40), max.edge = c(6, 9))
sigma0 = 0.05 ## field std.dev.
range0 = 20
kappa0 = sqrt(8)/range0
tau0 = 1/(sqrt(4*pi)*kappa0*sigma0)
spde = inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, 1), B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
                          theta.prior.prec = c(0.1, 1))


# Projection matrices
sp.mat = inla.spde.make.A(mesh, sp.xy)
quad.mat = inla.spde.make.A(mesh, quad.xy)
is.sp <- all.dat$Pres
is.quad <- 1 - all.dat$Pres


# Defining response
attach(all.dat)
stk.sp = inla.stack(data = list(y = 1, e = 0), A = list(sp.mat, 1),
                      tag='sp', effects = list(list(i = 1:mesh$n), data.frame(Intercept = 1,
                        FC = FC[is.sp], MNT = MNT[is.sp], MXT = MXT[is.sp], Rain = Rain[is.sp],
                        FC2 = FC[is.sp]^2, MNT2 = MNT[is.sp]^2, MXT2 = MXT[is.sp]^2,
                        Rain2 = Rain[is.sp]^2, FC.MNT = FC[is.sp]*MNT[is.sp],
                        FC.MXT = FC[is.sp]*MXT[is.sp], FC.Rain = FC[is.sp]*Rain[is.sp],
                        MNT.MXT = MNT[is.sp]*MXT[is.sp], MNT.Rain = MNT[is.sp]*Rain[is.sp],
                        MXT.Rain = MXT[is.sp]*Rain[is.sp], D.Main = sqrt(D.Main[is.sp]),
                        D.Urb = sqrt(D.Urb[is.sp]), D.Main2 = sqrt(D.Main[is.sp])^2,
                        D.Urb2 = sqrt(D.Urb[is.sp])^2,
                        D.Main.D.Urb = sqrt(D.Main[is.sp])*sqrt(D.Urb[is.sp]),
                        dist=0)))

stk.quad = inla.stack(data = list(y = 0, e = spat.res^2),
                        A = list(quad.mat, 1), tag='quad', effects = list(list(i = 1:mesh$n),
                          data.frame(Intercept = 1, FC = FC[is.quad], MNT = MNT[is.quad],
                          MXT = MXT[is.quad], Rain = Rain[is.quad], FC2 = FC[is.quad]^2,
                          MNT2 = MNT[is.quad]^2, MXT2 = MXT[is.quad]^2, Rain2 = Rain[is.quad]^2,
                          FC.MNT = FC[is.quad]*MNT[is.quad], FC.MXT = FC[is.quad]*MXT[is.quad],
                          FC.Rain = FC[is.quad]*Rain[is.quad], MNT.MXT = MNT[is.quad]*MXT[is.quad],
                          MNT.Rain = MNT[is.quad]*Rain[is.quad], MXT.Rain = MXT[is.quad]*Rain[is.quad],
                          D.Main = sqrt(D.Main[is.quad]), D.Urb = sqrt(D.Urb[is.quad]),
                          D.Main2 = sqrt(D.Main[is.quad])^2, D.Urb2 = sqrt(D.Urb[is.quad])^2,
                          D.Main.D.Urb = sqrt(D.Main[is.quad])*sqrt(D.Urb[is.quad]),
                          dist=1)))

stk.all = inla.stack(stk.sp, stk.quad)




# Fitting models
ft.inla = inla(y ~ 0 + Intercept + FC + MNT + MXT + Rain + FC2 + MNT2 + MXT2
               + Rain2 + FC.MNT + FC.MXT + FC.Rain + MNT.MXT + MNT.Rain + MXT.Rain + D.Main
               + D.Urb + D.Main2 + D.Urb2 + D.Main.D.Urb
               + f(inla.group(dist, n = 50, method = "quantile"), model = "rw1",
                   scale.model = TRUE) + f(i, model = spde), family = "poisson",
               data = inla.stack.data(stk.all),
               control.predictor = list(A = inla.stack.A(stk.all), compute = TRUE),
               E = inla.stack.data(stk.all)$e, control.compute = list(dic = TRUE),
               control.fixed = list(expand.factor.strategy="inla"))

# The fitted coefficients of ft.inla are presented in Table 3.
# We can project the mean and standard deviation of the latent field onto quadrature points
# using the inla.mesh.projector command as in Figure 4 of the main text:
proj.quad = inla.mesh.projector(mesh, quad.xy)
gf.mean = inla.mesh.project(proj.quad, ft.inla$summary.random$i$mean)
gf.sd = inla.mesh.project(proj.quad, ft.inla$summary.random$i$sd)

# We can produce a plot of the effect of the interaction covariate dist as follows:
plot(ft.inla$summary.random[[2]][, 1:2], type = 'l',
      xlab = 'dist (km)', ylab = 'f(dist)'); abline(h = 0, lty = 3)
for (i in c(4, 6)) lines(ft.inla$summary.random[[2]][,c(1, i)], lty = 2)


