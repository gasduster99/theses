rm(list=ls())

#
library(TTR)

#
#
#

#read.csv('ratioDat.csv')
#dat = read.csv('bh_1500_b0_3000_fm_1.5_M_050_300Dat.csv')
dat = read.csv('bRat_0.60_dRat_1.0_M_0.04_0.30DatSD.csv')
dat = dat[order(dat$M),]
head.x = dat$mean.x #dat$median.x
head.y = dat$mean.y #dat$median.y

#emaX = EMA(dat$head.x-dat$tail.x, n=4)
smaX5 = SMA(head.x-dat$tail.x, n=7)
fillX = c(
	(head.x-dat$tail.x)[1],
	SMA(head.x-dat$tail.x, n=2)[2],
        SMA(head.x-dat$tail.x, n=3)[3],
        #SMA(head.x-dat$tail.x, n=4)[4],
	na.omit(smaX5),	
        SMA(head.x-dat$tail.x, n=3)[length(head.x-dat$tail.x)-1],
	SMA(head.x-dat$tail.x, n=2)[length(head.x-dat$tail.x)-1],
       	(head.x-dat$tail.x)[length(head.x-dat$tail.x)]
)

#smaX  = SMA(dat$head.x-dat$tail.x, n=4)
pdf('bRat_0.60_dRat_1.0_M_0.04_0.30XBiasSD.pdf')
plot(dat$M, head.x-dat$tail.x,
	ylab="Bias(F*/M)",
	xlab="Natural Mortality"
)
lines(dat$M, fillX, lwd=3)
#lines(dat$M, emaX)
#lines(dat$M, smaX5, lwd=3)
#lines(dat$M[1:5], fillX, lwd=3)
#lines(dat$M, smaX, col='blue')
##lines(dat$M, vmaX, col='green')
#lines(dat$M, wmaX, col='green')
##
dev.off()
#
pdf('bRat_0.60_dRat_1.0_M_0.04_0.30YBiasSD.pdf')
#dev.new()
#emaY = EMA(dat$head.y-dat$tail.y, n=4)
smaY5 = SMA(head.y-dat$tail.y, n=7)
#fillY = c(
#        (head.y-dat$tail.y)[1],
#        SMA(head.y-dat$tail.y, n=2)[2],
#        #SMA(head.y-dat$tail.y, n=3)[3],
#        #SMA(head.y-dat$tail.y, n=4)[4],
#        na.omit(smaY5)[1],
#	SMA(head.y-dat$tail.y, n=2)[length(head.y-dat$tail.y)-1],
#	(head.y-dat$tail.y)[length(head.y-dat$tail.y)]
#)
#smaY  = SMA(dat$head.y-dat$tail.y, n=4)
fillY = c( 
	(head.y-dat$tail.y)[1],
	SMA(head.y-dat$tail.y, n=2)[2],
	SMA(head.y-dat$tail.y, n=3)[3],
	na.omit(smaY5),
	SMA(head.y-dat$tail.y, n=3)[length(head.y-dat$tail.y)-2],
	SMA(head.y-dat$tail.y, n=2)[length(head.y-dat$tail.y)-1],
       (head.y-dat$tail.y)[length(head.y-dat$tail.y)]
)


plot(dat$M, head.y-dat$tail.y,
	ylab=expression("Bias("*B[h]/B[0]*")"),
	xlab="Natural Mortality"
)
#lines(dat$M, emaY)
smaY5[1] = (head.y-dat$tail.y)[1]
lines(dat$M, fillY, lwd=3)
#lines(dat$M, smaY5, lwd=3)
#lines(dat$M[1:5], fillY, lwd=3)
#lines(dat$M, smaY, col='blue')
##lines(dat$M, vmaX, col='green')
#lines(dat$M, wmaY, col='green')
dev.off()
