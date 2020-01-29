rm(list=ls())

#
library(TTR)

#
#
#

#
files = Sys.glob('dRat_*_M_0.04_0.30Dat.csv')
#
Ms = list()
yFills = list()
xFills = list()
#
for(file in files){
	dat = read.csv(file)
	dat = dat[order(dat$M),]
	head.x = dat$median.x #dat$mean.x 
	head.y = dat$median.y #dat$mean.y 
	#
	smaX5 = SMA(head.x-dat$tail.x, n=5)
	fillX = c(
		(head.x-dat$tail.x)[1],
		SMA(head.x-dat$tail.x, n=2)[2],
	        SMA(head.x-dat$tail.x, n=3)[3],
	        SMA(head.x-dat$tail.x, n=4)[4],
		na.omit(smaX5) #na.omit(smaX5)[1]
	)
	#
	smaY5 = SMA(head.y-dat$tail.y, n=5)
	fillY = c(
	        (head.y-dat$tail.y)[1],
	        SMA(head.y-dat$tail.y, n=2)[2],
	        SMA(head.y-dat$tail.y, n=3)[3],
	        SMA(head.y-dat$tail.y, n=4)[4],
	        na.omit(smaY5)#na.omit(smaY5)[1]
	)
	#
	Ms[[file]] = dat$M
	xFills[[file]] = fillX
	yFills[[file]] = fillY
}
#
cols = rep(brewer.pal(11, 'Spectral')[-6], 2)
#
dev.new()
plot(Ms[[1]], yFills[[1]], 'l', 
	ylim=c(do.call(min, yFills), do.call(max, yFills)),
	xlim=c(do.call(min, Ms), do.call(max, Ms)),
	col=cols[1]
)
i = 2
for(f in files[-1]){ lines(Ms[[f]],yFills[[f]], col=cols[i]); i=i+1; }
#
dev.new()
plot(Ms[[1]], xFills[[1]], 'l', 
	ylim=c(do.call(min, xFills), do.call(max, xFills)),
	xlim=c(do.call(min, Ms), do.call(max, Ms)),
	col=cols[1]
)
i = 2
for(f in files[-1]){ lines(Ms[[f]],xFills[[f]], col=cols[i]); i=i+1; }


#files = c(
#	'bh_1200_b0_3000_5M_040_300Dat.csv', 
#	'bh_1500_b0_3000_5M_040_300Dat.csv', 
#	'bh_1800_b0_3000_5M_040_300Dat.csv', 
#	'bh_2700_b0_3000_5M_040_300Dat.csv'
#)
##files = c(
##	'bh_0900_b0_3000_M_040_300Dat.csv', 
##	'bh_1200_b0_3000_M_040_300Dat.csv', 
##	'bh_1500_b0_3000_M_040_300Dat.csv', 
##	'bh_1800_b0_3000_M_040_300Dat.csv',
##)
#Ms = list(NA, NA, NA, NA)
#names(Ms) = files
#yFills = list(NA, NA, NA, NA)
#names(yFills) = files
#xFills = list(NA, NA, NA, NA)
#names(xFills) = files

	##smaX  = SMA(dat$head.x-dat$tail.x, n=4)
	#pdf('bh_2400_b0_3000_M_040_300XBias.pdf')
	#plot(dat$M, head.x-dat$tail.x,
	#	ylab="Bias(F*/M)",
	#	xlab="Natural Mortality"
	#)
	##lines(dat$M, emaX)
	#lines(dat$M, smaX5, lwd=3)
	#lines(dat$M[1:5], fillX, lwd=3)
	##lines(dat$M, smaX, col='blue')
	###lines(dat$M, vmaX, col='green')
	##lines(dat$M, wmaX, col='green')
	##
	#dev.off()
	##
	#pdf('bh_2400_b0_3000_M_040_300YBias.pdf')
	##dev.new()
	##emaY = EMA(dat$head.y-dat$tail.y, n=4)
	##smaY  = SMA(dat$head.y-dat$tail.y, n=4)
	#plot(dat$M, head.y-dat$tail.y,
	#	ylab=expression("Bias("*B[h]/B[0]*")"),
	#	xlab="Natural Mortality"
	#)
	##lines(dat$M, emaY)
	#smaY5[1] = (head.y-dat$tail.y)[1]
	#lines(dat$M, smaY5, lwd=3)
	#lines(dat$M[1:5], fillY, lwd=3)
	##lines(dat$M, smaY, col='blue')
	###lines(dat$M, vmaX, col='green')
	##lines(dat$M, wmaY, col='green')
	#dev.off()

