rm(list=ls())

library(tgp)
library(akima)

name = "Lock240000"; uu=50; lambda=0.25

#name = "RosegpPic"; uu=30; tMin=0; lambda=0.2

#name = "RastHard"; uu=40; tMin=0; lambda=0.3

#name = "EasomMed"; uu=20; tMin=-1; lambda=0.1


##name = "RoseEasyEasy"; uu=30; tMin=0; lambda=0.2
##name = "RastHard2"; uu=30; tMin=0; lambda=0.3


load(sprintf("seaSave%s.RData", name))
W = uu
tMin = min(ZMax = seaSave[[1]][[M]])

#
out = init[[6]]
figPlace = "/home/nick/Documents/school/ucscGrad/thesis/herbie/writting/figures/"

##
##pdf(sprintf("%sgpMean%s.pdf", figPlace, name))
##
#plot(out$obj, 
#	layout="surf", 
#	phi=     60,    #20,  #40. 
#	theta=   90    #-110  #-40
#)
##dev.off()

#dev.new()
pdf(sprintf("%sgpMean%s.pdf", figPlace, name))
tmp = interp(init[[6]]$obj$XX$x1, init[[6]]$obj$XX$x2, init[[6]]$obj$ZZ.mean)
#
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
filled.contour( tmp[[1]], tmp[[2]], tmp[[3]],#init[[6]]$obj$XX$x1, init[[6]]$obj$XX$x2, init[[6]]$obj$ZZ.mean,
                        color.palette = jet.colors,
                        xlab=expression(x[1]), #expression(Q[A[1]]),
                        ylab=expression(x[2]), #expression(Q[A[2]]),
                        main='z mean',
			xlim=c(0, 20000),
			ylim=c(0, 20000)
                        #zlim=c(0, 3000)
#                       plot.axes={ points(c(pi), c(pi), col="white", pch=20); axis(1); axis(2) }
#
)
dev.off()

#dev.new()
#plot(out$obj, 
#	as = "improv",
#	layout="as"
#)


