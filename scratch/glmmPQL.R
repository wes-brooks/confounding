library(nlme)
library(MASS)


cs1Exp <- corExp(0.03, form=~x+y, fixed=TRUE)
cs1Exp <- Initialize(cs1Exp, loc.tot)

cs1

D.tot <- Matrix(0, n.tot, n.tot)
for (i in 1:(n.tot-1)) {indx <- dd.tot$delsgs$ind1[dd.tot$delsgs$ind2==i]; if (length(indx) > 0) D.tot[i,indx] <- 1}

ft.glmmPQL <- glmmPQL(resp~XX, data=dat, family='poisson', weights=wt, random=~1|id, correlation=cs1Exp)
