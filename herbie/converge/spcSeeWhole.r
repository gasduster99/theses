rm(list=ls())

##load("seaSaveRastHard.RData")
##load("seaSaveRoseSkew.RData")



#load("seaSaveRastEasy.RData")
#load("seaSaveRoseBig.RData")
#load("seaSaveRoseEasy.RData")
load("seaSaveRoseEasyEasy.RData")

 
sig = 1#2#3#2.5
W = 10 

#
wSamples = c()
#
dev.new(width=14, height=14/3)
#
it = as.numeric(readline("Start Where? "))

while(T){
	ZMax = seaSave[[1]][[it]]
	samples = seaSave[[2]][[it]]	
	maxes = seaSave[[3]][[it]]
	mcmcSize = seaSave[[4]][[it]]
	#
	sEnd = length(samples)
	back = W*mcmcSize
	lEnd = max(0, sEnd-back)
	wSamples = samples[seq(lEnd, sEnd)] 
	#
        wm = mean(wSamples)
        wv = var(wSamples)
        wmu = log( (wm^2)/(wv+wm^2)^0.5 )
        ws2 = log( 1 + (wv/(wm^2)) )
        wUpper3 = exp( wmu + (sig*(ws2^0.5)) )
        wLower3 = exp( wmu - (sig*(ws2^0.5)) )
        #
        samLower = quantile(wSamples, probs=pnorm(-sig))
        samUpper = quantile(wSamples, probs=pnorm(sig))
	prob = pnorm(sig) - pnorm(-sig)
	#
	#OUTPUT
        #ding(1)
        layout( matrix(c(1, 2, 3), nrow=1, ncol=3) )
        #
        #Figure 1
        plot( seq(1, it), ZMax,
                "l",
                xlab="Sweeps",
                ylab="Value",
                main="Best Z Value"
        )
        #
	it=it-1
	#
        #Figure 2
        L = length(seq(1, it))
        left = max((L-(W-1)), 1)
	wOut = maxes[1:it]<wLower3 | maxes[1:it]>wUpper3
	samOut = maxes[1:it]<samLower | maxes[1:it]>samUpper
        plot( seq(1, it), maxes[1:it],
                main=sprintf('Sweeps: %d\n', it),
                ylab='Max Improvement\n',
                col=c(rep("black", L-1), "green"),
                ylim=c(0, max(wUpper3, maxes[1:it], samUpper))  #0.007)
        )
	points( seq(1, it)[wOut], maxes[1:it][wOut], 
		pch="O",
		col="red",
		cex=1.3
	)
	points( seq(1, it)[samOut], maxes[1:it][samOut], 
		pch="O",
		col="blue",
		cex=2
	)
        #Sample Intervals
        segments(    1, samUpper, left, samUpper, col="blue", lty=2)
        segments(    1, samLower, left, samLower, col="blue", lty=2)
        segments( left, samUpper,        it, samUpper, col="blue")
        segments( left, samLower,        it, samLower, col="blue")
        #log(Normal Intervals)
        segments(    1, wUpper3, left, wUpper3, col="red", lty=2)
        segments(    1, wLower3, left, wLower3, col="red", lty=2)
        segments( left, wUpper3,   it, wUpper3, col="red")
        segments( left, wLower3,   it, wLower3, col="red")
        #BarBar
        segments( left, wm, it, wm, col="black")
        #
        #Figure 3
        hist( wSamples, add=F, freq=F,
                main="Window Samples"
        )
        abline( v=wUpper3, col="red")
        abline( v=wLower3, col="red")
        abline( v=samUpper, col="blue")
        abline( v=samLower, col="blue")
        abline( v=maxes[it], col='green' )
        abline( v=wm, col='black' )
        legend("topright",
                legend=c(sprintf("%2.3f%s (%d sigma)\nLog-N Interval\n", prob*100, "%", sig), 
			 sprintf("%2.3f%s (%d sigma)\nSample Interval\n", prob*100, "%", sig), 
			 "x-bar", 
			 "x-bar-bar"),
                col=c("red", "blue", "green", "black"),
                lty=1
        )

	move = readline("")

	if( move=="d" | move=="w" ){
		it = it+2
	}else if( move=="a" | move=="s" ){
		it = it
	}else{
		it = as.numeric(move) + 1
	}
	
}

