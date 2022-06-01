rm(list=ls())

#
library(qcc)
#
source("spcSeeEWMANormFunc.r")

#
#

#name = "RoseEasyEasy"; uu=30; tMin=0;  	#it=74          #lambda=0.2
name = "RastHard"; uu=59; tMin=0;	#it=138
#uu=40; tMin=0;	#it=96          #lambda=0.3

#
#

load(sprintf("seaSave%s.RData", name))
W = uu
it = 138
h = 0.01

#
#

it0 = W+1
it = it0
while(T){
	#
	ewmaOut = findEWMA(it, h)
	ewmaSS = findSS(it, h)
	#
	pdf(sprintf("ewmaConvChart%s.pdf", name))
	ewmaPlot(ewmaOut, it)
	legend( "topright", c(expression(Y[i]), expression(Z[i]*" Violation"), expression(Z[i]*" Control"), "Control Limit", NA,"\nControl Window\nBoundary"), 
	               pch=c(3 , 1, 20, NA, NA, NA),
	               lty=c(NA, NA, NA, 2, 3, NA),
	               col=c("black", "red", "black", "black", "black", NA)
	       )
	dev.off()
	#
	pdf(sprintf("ss%s.pdf", name))
	plot(ewmaSS[[1]], ewmaSS[[2]], 'l',
		xlab = expression(lambda),
                ylab = expression('Forcasting S'[lambda]),
		main = expression('Finding the Optimal '*hat(lambda)) 
	)
	abline(v=ewmaOut[[10]], col='red', lty=2)
	dev.off()
	#
	pdf(sprintf("bestZ%s.pdf", name))
	ZMax = ewmaOut[[11]]
        plot( seq(1, it), ZMax,
                "l",
                xlab="Iteration Number",
                ylab="Objective Value",
                main="Best Z Value",
                ylim=c(tMin, max(ZMax))
        )
        abline( h=tMin, lty=3 )
        legend("left", c("\nCurrent Best\n Z Value", "\nTheoretical\n Minimum", NA),
                lty=c(1, 3, NA)
        )
	dev.off()
	#
	it = wasd(name)
}











#	#if( it!=it0 ){ dev.off() }
#	#
#	ZMax = seaSave[[1]][[it]]
#	samples = seaSave[[2]][[it]]	
#	#maxes = seaSave[[3]][[it]]
#	mcmcSize = seaSave[[4]][[it]]	
#	#
#	sEnd = length(samples)
#        back = W*mcmcSize
#        lEnd = max(1, sEnd-back)
#        wSamples = samples[seq(lEnd, sEnd)]
#	mSamples = slice(samples, mcmcSize) 
#	#
#	wOut = LNToN(wSamples)
#	wm = wOut[[1]][1]
#	wv = wOut[[1]][2]
#	wmu = wOut[[2]][1]
#	ws2 = wOut[[2]][2]
#	#
#	mOut = lapply(mSamples, LNToN)
#	mm = getVect(mOut, "LN", "M") 
#	mv = getVect(mOut, "LN", "V")
#	mmu = getVect(mOut, "N", "M")
#	ms2 = getVect(mOut, "N", "V")	
#	#
#	#Figure 1
#	pdf(sprintf("bestZ%s.pdf", name))
#	plot( seq(1, it), ZMax,
#                "l",
#                xlab="Iteration Number",
#                ylab="Objective Value",
#                main="Best Z Value", 
#		ylim=c(tMin, max(ZMax))
#        )
#	abline( h=tMin, lty=3 )
#	legend("left", c("\nCurrent Best\n Z Value", "\nTheoretical\n Minimum", NA),
#		lty=c(1, 3, NA)
#	)
#	dev.off()
#	#
#	ewmaOut = ewma( rev(mmu)[1:W], lambda=lambda, newdata=rev(mmu)[W:length(mmu)], plot=F )	
#	testOut = rev(test1())
#	fudge = sample(rev(sort(ewmaOut$statistics))[1:3], 1)
#	min2 = min(mmu, ewmaOut$y, ewmaOut$limits)
#	max2 = max(mmu, ewmaOut$y, ewmaOut$limits)
#	#
#	#Figure 2
#	#dev.new()
#	pdf(sprintf("ewmaConvChart%s.pdf", name))
#	plot( seq(1, it), rev(c(ewmaOut$statistics, ewmaOut$newstats)),
#		pch=3, 
#		ylim=c(min2, max2),
#		xlab="Iteration Number",
#		ylab="Summary Statistics",
#		main="EWMA Convergence Chart"
#	)
#	lines( ewmaOut$x, rev(ewmaOut$y) )
#	lines( seq(1, it), rev(ewmaOut$limits[,1]), lty=2 )
#	lines( seq(1, it), rev(ewmaOut$limits[,2]), lty=2 )
#	points( ewmaOut$x, rev(ewmaOut$y), pch=20 ) 
#	points( ewmaOut$x[testOut], rev(ewmaOut$y)[testOut], pch=19, col="red" )		
#	abline( v=(it+0.5)-W, lty=3 )
#	abline( h=ewmaOut$center )
#	legend( "topright", c(expression(Y[i]), expression(Z[i]*" Violation"), expression(Z[i]*" Control"), "Control Limit", NA,"\nControl Window\nBoundary"), 
#		pch=c(3 , 19, 20, NA, NA, NA),
#		lty=c(NA, NA, NA, 2, 3, NA),
#		col=c("black", "red", "black", "black", "black", NA)
#	)
#	dev.off()
#	
#	#
#	#it = it-1
#	#	
#	it = wasd(name)	
#}































