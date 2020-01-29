rm(list=ls())

#
#library(TTR)
library(RColorBrewer)
#
source('dotCurve.r')

#
#FUNCTIONS
#

#
centerSMA = function(x, n){
	#x: values to be smoothed
	#n: an odd sized(n>=3) size of smoothing window
	
	#
	library(TTR)
	
	#
	l = length(x)
	m = 2:(floor(n/2))
	ml = length(m)
	#
	out = matrix(NA, nrow=length(x), ncol=1)
	out[1] = x[1]	
	for(i in m){ out[i]=SMA(x, n=i)[i] }
	#first+m+next
	out[(1+ml+1):(l-1-ml)] = na.omit(SMA(x, n=n))
	for(i in rev(m)-1){ out[l-i]=SMA(x, n=i)[l-i] } 
	out[l] = x[l]
	#
	return(out)
}

#
#MAIN
#

#
files = Sys.glob('bRat_*_dRat_*_M_0.04_0.30DatSD.csv')
#
Ds = list()
Cs = list()
Ms = list()#NA, NA, NA, NA)
#names(Ms) = files
yFills = list()#NA, NA, NA, NA)
#names(yFills) = files
xFills = list()#NA, NA, NA, NA)
#names(xFills) = files
for(f in files){
	#
	c = as.numeric(strsplit(f, '_')[[1]][c(2,4)])
	#
	dat = read.csv(f)
	dat = dat[order(dat$M),]
        head.x = dat$median.x  
        head.y = dat$median.y
	#only stuff that ran at least a little
	if( dim(dat)[1]>0 ){
		#
		x = head.x-dat$tail.x
		xHat = centerSMA(x, n=5)
		#
		y = head.y-dat$tail.y
                yHat = centerSMA(y, n=5)
		#
        	Ms[[f]] = dat$M
        	xFills[[f]] = xHat
        	yFills[[f]] = yHat		
		#
		Cs[[f]] = c
		Ds[[f]] = missSpec(c[1], c[2], c(0, 5))
	}
}

#
#PLOT
#

#
fs = names(Ms)
#cols = c('black', rep(brewer.pal(9, 'Set1'), 3))#rep(brewer.pal(length(Ms$bRat_0.60_dRat_1.0_M_0.04_0.30DatSD.csv), 'Spectral'), 3))
#
plot(Ds[[1]], yFills[[1]][1], 
	xlim=c(0, 2.5), 
	ylim=c(-0.6, 0)
	#cols=cols[1],
)
for(f in fs){
	for(j in 1:length(Ms[[f]])){
		#	
		points(Ds[[f]], yFills[[f]][j] 
			#cols=cols[j]
		)
	}
}
#
dev.new()
plot(-1, -1, #Ds[[1]], xFills[[1]][1], 
	xlim=c(0, 2.5), 
	ylim=c(-0.6, 30)
	#cols=cols[1],
)
for(j in 1:24){#length(Ms[[f]])){
	ys = c()
	ds = c()
	for(f in fs){		
		points(Ds[[f]], xFills[[f]][j]
			#cols=cols[j]
		)
                ys = c(ys, xFills[[f]][j])
		ds = c(ds, Ds[[f]])
	}
	print('')
	print(length(ys))
	print(length(ds))
}







#plot(Ds[[fs[1]]], xFills[[fs[1]]][1], xlim=c(0, 2.5), ylim=c(-0.5, 30))
#for(f in fs[-1]){ points(Ds[[f]], xFills[[f]][1]) }










#
#BONE YARD
#

	##
	#smaX5 = SMA(head.x-dat$tail.x, n=7)
	#fillX = c(
        #(head.x-dat$tail.x)[1],
        #SMA(head.x-dat$tail.x, n=2)[2],
        #SMA(head.x-dat$tail.x, n=3)[3],
        ##SMA(head.x-dat$tail.x, n=4)[4],
        #na.omit(smaX5),
        #SMA(head.x-dat$tail.x, n=3)[length(head.x-dat$tail.x)-1],
        #SMA(head.x-dat$tail.x, n=2)[length(head.x-dat$tail.x)-1],
        #(head.x-dat$tail.x)[length(head.x-dat$tail.x)]

##
#x = 1:20
#y = rnorm(20, sin(x), 1)
##
#yHat = centerSMA(y, n=7)
##
#plot(x, y)
#lines(x, yHat, lwd=3)


