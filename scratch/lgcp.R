polygon <- data.frame(x=c(0,1,0,1), y=c(0,0,1,1))
sd <- ppp(loc.cox$x, loc.cox$y, poly=polygon)

minimum.contrast(sd, model='exponential', method='g', intens=density(sd), transform=log)
chooseCellwidth(sd, cwinit=0.03)
Cellwidth <- 0.03

n.side <- 50
quad = data.frame(x=rep(seq(0,1,n.side + 2)[], times=n.side), y=rep(seq(0,1,n.side), each=n.side))
covar <- SpatialPixelsDataFrame(cbind(quad$x, quad$y), quad)
polyolay <- getpolyol(data=sd, pixelcovariates=covar, cellwidth=Cellwidth)
