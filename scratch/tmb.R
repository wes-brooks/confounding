X2 <- as.matrix(cbind(1, (XX-mean(XX))/sd(XX)))

#Q <- as.matrix(t(SRE.ortho) %*% Precision %*% SRE.ortho)
Q <- as.matrix(Precision)
logdetQ <- eigen(Q)$values %>% log %>% sum
SRE.local <- as.matrix(SRE.ortho)
n.random <- ncol(SRE.local)

# # 
library(TMB)
setwd("src")
compile("spatial_reduced.cpp")
dyn.load(dynlib("spatial_reduced"))
dyn.load(dynlib("spatial_penalized"))
obj <- MakeADFun(data=list(
                   y=dat$resp,
                   X=X2,
                   Q=Q,
                   logdetQ=logdetQ,
                   SRE=SRE.local,
                   wt=p.wt,
                   ltau = 0.0
                 ),
                 parameters=list(
                   #b=c(log(sum(p.wt*dat$resp / sum(p.wt))),0),
                   b = c(0,0),
                   #ltau=0.0,
                   u = rep(0, n.random)),
                 random=c("u"),
                 DLL="spatial_reduced")
obj$control <- list(trace=2, reltol=1e-2, maxit=100)
time.fit <- system.time(opt<-optim(obj$par, obj$fn, obj$gr, method='CG', control=obj$control))
rep <- sdreport(obj)
rep

sum(rep$par.random)
sum(rep$par.random * XX)

sourceCpp("src/ldl.cpp")
logdetQ <- ldl(Q2) %>% `[`(-n.tot) %>% log %>% sum
u = rep$par.random
ltau = opt$par[['ltau']]
eta = opt$par[['a']] + opt$par[['b']] * XX + u
lam = exp(eta)
sum(dpois(dat$resp, lam, log=TRUE) %>% `*`(p.wt)) - logdetQ + exp(ltau)/2 * (t(u) %*% Q2 %*% u) + n.tot/2.0 * log(2*pi - ltau)


obj <- MakeADFun(data=list(y=dat$resp, X=X2, Q=Q, L=LIE$vectors[,1:5], wt=p.wt, ridge=0.5),
                 parameters=list(
                   b=c(log(sum(p.wt*dat$resp / sum(p.wt))), 0),
                   ltau=0.0,
                   u = rep(0, 5)
                 ),
                 random=c("u"),
                 DLL="spatial_reduced")
