rm(list=ls())

library(tgp)
library(akima)

colPersp = function(x,y,z,pal,nb.col,...,xlg=TRUE,ylg=TRUE){
        colnames(z) <- y
        rownames(z) <- x

        nrz <- nrow(z)
        ncz <- ncol(z)

        color <- pal(nb.col)
        zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
        facetcol <- cut(zfacet, nb.col)
        par(xlog=xlg,ylog=ylg)
        persp(
                as.numeric(rownames(z)),
                as.numeric(colnames(z)),
                as.matrix(z),
                col=color[facetcol],
                ...
                )
}

name = "RosegpPic"; uu=30; tMin=0; lambda=0.2
#name = "RoseEasyEasy"; uu=30; tMin=0; lambda=0.2

#name = "RastHard"; uu=40; tMin=0; lambda=0.3
##name = "RastHard2"; uu=30; tMin=0; lambda=0.3

#name = "EasomMed"; uu=20; tMin=-1; lambda=0.1

load(sprintf("seaSave%s.RData", name))
W = uu

#
out = init[[6]]
figPlace = "/home/nick/Documents/school/ucscGrad/thesis/herbie/writting/figures/"

#
#pdf(sprintf("%sgpMean%s.pdf", figPlace, name))
tmp = interp(init[[6]]$obj$XX$x1, init[[6]]$obj$XX$x2, init[[6]]$obj$ZZ.mean)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colPersp( tmp[[1]], tmp[[2]], tmp[[3]],
        jet.colors,
        200,
        phi=45,
        theta=-40,
        border=NA,
        box=T,
        ticktype="detailed",
        xlab="x1",
        ylab="x2",
        #zlab='Rosenbrock'
        zlab="z mean"
        #main='Rosenbrock'
)

##
#plot(out$obj, 
#	layout="surf", 
#	phi=     40, 
#	theta=   -40
#)
#dev.off()

#
#dev.new()
#plot(out$obj, 
#	as = "improv",
#	layout="as"
#)


