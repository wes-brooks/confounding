library(ggplot2)
library(dplyr)
library(plyr)

tiles <- tile.list(dd.tot)
tiles <- lapply(tiles, function(p) as.data.frame(p[c("x","y")]))


for (i in 1:length(tiles)) colnames(tiles[[i]] < c('x', 'y'))
tiles %>% sapply(function(x) Polygon(x)) -> polylist
(1:length(polylist)) %>% sapply(function(i) Polygons(list(polylist[[i]]), ID=i)) %>% SpatialPolygons -> spoly
SpatialPolygonsDataFrame(spoly, data=dat) -> spdf

dd <- cbind(dat, b, rando)

colnames(S) <- paste("SRE", 1:ncol(S), sep='')
dat2 <- cbind(dat, as.matrix(S))

colnames(Seff) <- paste("SRE", 1:ncol(Seff), sep='')
dat2 <- cbind(dat, as.matrix(Seff))

SRE.ortho <- P.c %*% SRE
dat3 <- cbind(dat, as.matrix(SRE.ortho))

SRE.par <- P %*% SRE
dat4 <- cbind(dat, as.matrix(SRE.par))

spdf.flat <- fortify(spdf)
spdf.plot <- join(spdf.flat, dd, by='id')
ggplot(spdf.plot) + aes(x=long, y=lat, group=group, fill=b) + geom_polygon() + scale_fill_gradient2() + geom_point(data=loc.cox, aes(group=NULL, fill=NULL))


spdf.flat.2 <- fortify(spdf)
spdf.plot.2 <- join(spdf.flat.2, dat2, by='id')
ggplot(spdf.plot.2) + aes(x=long, y=lat, group=group, fill=SRE52) + geom_polygon() + scale_fill_gradient2()


spdf.flat.ortho <- fortify(spdf)
spdf.flat.ortho <- join(spdf.flat.ortho, dat3, by='id')
ggplot(spdf.flat.ortho) + aes(x=long, y=lat, group=group, fill=SRE2) + geom_polygon() + scale_fill_gradient2()


spdf.flat.parallel <- fortify(spdf)
spdf.flat.parallel <- join(spdf.flat.parallel, dat4, by='id')
ggplot(spdf.flat.parallel) + aes(x=long, y=lat, group=group, fill=SRE2) + geom_polygon() + scale_fill_gradient2()
