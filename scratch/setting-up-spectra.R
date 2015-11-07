X2 <- as.matrix(cbind(1,XX))

P.1 <- X2 %*% solve(t(X2) %*% X2)
P.2 <- t(X2) %*% SRE

SRE.ortho <- SRE - P.1 %*% P.2

P <- X2 %*% solve(t(X2) %*% X2) %*% t(X2)
P.c <- diag(n.tot) - P


IA <- P.c %*% off.diag %*% P.c
LIE <- eigen(IA)


L <- eigen(P.c)$vectors[,1:(n.tot-1)]

QQ <- as.matrix(t(L) %*% Q %*% L)
