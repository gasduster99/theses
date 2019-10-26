rm(list=ls())

library(qcc)

LNToN = function(wm, wv){
        #wm = mean(LNSamples)
        #wv = var(LNSamples)
        #associated normal stats
        wmu = log( (wm^2)/(wv+wm^2)^0.5 )
        ws2 = log( 1 + (wv/(wm^2)) )
        #
        #           LN Stats    N Stats
        out = list( c(wm, wv), c(wmu, ws2) )
}

test1 = function(int){
        out = ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1] 
        #
        return(out)
}


#
#tune that shit
fFactor = 1000
a = 0
c = 69
bOne = 10
bTwo = 55
#
domain = seq(a, c)
out = matrix(NA, nrow=length(domain), ncol=1)
#
#create the data
i=1
for(x in domain){
	if(x < bOne){
		out[i] = (10+rnorm(1, 0, 5))/fFactor		
	}else if(bOne <= x & x < bTwo){
		out[i] = (4+rnorm(1, 0, 1))/fFactor
	}else{
		out[i] = (0+rbeta(1, 2, 5)*175000/((4*i)^2))/fFactor
	}
	i=i+1
}
#
#EWMA
mmu = out
W = 20
it=length(out)+1
#
ewmaOut = ewma( rev(mmu)[1:W], newdata=rev(mmu)[W:length(mmu)], plot=F )        
testOut = rev(test1())
min2 = min(mmu, ewmaOut$y, ewmaOut$limits)
max2 = max(mmu, ewmaOut$y, ewmaOut$limits)

#
#example EWMA Convergence Chart
pdf("exampleEWMA.pdf")
plot( seq(1, it), rev(c(ewmaOut$statistics, ewmaOut$newstats)),
        pch=3,
        ylim=c(min2, max2),
	xlim=c(0, c+4),
        xlab="Iteration Number",
        ylab="Summary Statistics",
        main="EWMA Convergence Chart"
)
lines( ewmaOut$x, rev(ewmaOut$y) )
lines( seq(1, it), rev(ewmaOut$limits[,1]), lty=2 )
lines( seq(1, it), rev(ewmaOut$limits[,2]), lty=2 )
points( ewmaOut$x, rev(ewmaOut$y), pch=20 )
points( ewmaOut$x[testOut], rev(ewmaOut$y)[testOut], pch=19, col="red" )
abline( v=(it+0.5)-W, lty=3 )
abline( h=ewmaOut$center )
legend("topright", c(expression(Y[i]), expression(Z[i]*" Violation"), expression(Z[i]*" Control"), "Control Limit", NA,"\nControl Window\nBoundary"),
        pch=c(3 , 19, 20, NA, NA, NA),
        lty=c(NA, NA, NA, 2, 3, NA),
        col=c("black", "red", "black", "black", "black", NA)
)
dev.off()
#
#Example Progression
pdf("exampleEI.pdf")
plot(domain, out, 
	xlab="Iteration", 
	ylab="E[ I(x) ]", 
	main="Example Convergence"
)
rect(a-bOne, 0-bOne, bOne, max(out)+bOne, 
	col="lightgrey",
	border=NA
)
rect(bOne, 0-bOne, bTwo, max(out)+bOne, 
	col="grey",
	border=NA
)
rect(bTwo, 0-bOne, c+bOne, max(out)+bOne, 
	col="darkgrey",
	border=NA	
)
points(domain, out, 
	xlab="Iteration", 
	ylab="E[ I(x) ]", 
	main="Example Convergence", 
)
dev.off()
#
load("./converge/seaSaveRoseEasy.RData")
it=110#100#80#70#50#20
#
#Example Improvement Histogram
pdf("exampleIHist.pdf")
hOut = hist( seaSave[[2]][[it]], 
	main="Example Improvement Criteria Histogram", 
	xlab="Improvement Criteria"
)
#points(mean(seaSave[[2]][[it]]), hOut$counts[hOut$breaks<mean(seaSave[[2]][[20]])],
#	col="red"
#)
segments(mean(seaSave[[2]][[it]]), 0, mean(seaSave[[2]][[it]]), hOut$counts[hOut$breaks<mean(seaSave[[2]][[it]])],
	col="blue"
)
legend("topright", legend=c("Expected Value"), col=c("blue"), lty=c(1))
dev.off()
#
#Example log-Improvement Histogram
pdf("exampleLogIHist.pdf")
hOut = hist( log(seaSave[[2]][[it]][seaSave[[2]][[it]]!=0]), 
        main="Example Log-Improvement Criteria Histogram",
        xlab="Log-Improvement Criteria"
)
segments(mean(log(seaSave[[2]][[it]][seaSave[[2]][[it]]!=0])), 0, mean(log(seaSave[[2]][[it]][seaSave[[2]][[it]]!=0])), hOut$counts[hOut$breaks<mean(log(seaSave[[2]][[it]][seaSave[[2]][[it]]!=0]))],
        col="blue"
)
legend("topleft", legend=c("Expected Value"), col=c("blue"), lty=c(1))
dev.off()

