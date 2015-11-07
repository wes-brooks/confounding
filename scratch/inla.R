# Running simulations through INLA

library(INLA)

hull <- inla.nonconvex.hull(as.matrix(loc.tot), resolution=100)
mesh <- inla.mesh.2d(boundary=hull, cutoff=30/sqrt(86000), max.edge=c(6,9)/sqrt(86000))
sigma0 <- 0.05
range0 <- 20 / sqrt(86000)
kappa0 <- sqrt(8) / range0
tau0 <- 1 / sqrt(4*pi) / kappa0 / sigma0

spde <- inla.spde2.matern(mesh, B.tau=cbind(log(tau0), -1, 1),
                          B.kappa=cbind(log(kappa0), 0, -1),
                          theta.prior.mean=c(0,0),
                          theta.prior.prec=c(0.1, 1))

sp.mat <- inla.spde.make.A(mesh, as.matrix(loc.cox))
quad.mat <- inla.spde.make.A(mesh, as.matrix(loc.q))

stack.cox <- inla.stack(data=list(y=1, e=0),
                        A=list(sp.mat, 1),
                        tag='sp',
                        effects=list(list(i=1:mesh$n),
                                    data.frame(Intercept=1,
                                               XX=X.cox)
                                    )
                        )

stack.q <- inla.stack(data=list(y=0, e=1/105),
                        A=list(quad.mat, 1),
                        tag='sp',
                        effects=list(list(i=1:mesh$n),
                                    data.frame(Intercept=1,
                                               XX=X.q)
                                    )
                        )

stack.all <- inla.stack(stack.cox, stack.q)

ft.inla <- inla(y ~ 0 + Intercept + XX + f(i, model=spde),
                family='poisson',
                data=inla.stack.data(stack.all),
                control.predictor=list(A=inla.stack.A(stack.all), compute=TRUE),
                E=inla.stack.data(stack.all)$e,
                control.compute=list(dic=TRUE),
                control.fixed=list(expand.factor.strategy='inla'))


proj.quad <- inla.mesh.projector(mesh, as.matrix(loc.q))
proj.all <- inla.mesh.projector(mesh, as.matrix(loc.tot))

gf.mean <- inla.mesh.project(proj.all, ft.inla$summary.random$i$mean)
gf.sd <- inla.mesh.project(proj.quad, ft.inla$summary.random$i$sd)

