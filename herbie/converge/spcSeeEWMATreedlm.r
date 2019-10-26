rm(list=ls())

library(tgp)

#load("seaSaveRoseEasyEasy.RData")
#load("seaSaveRoseEasy.RData")
#load("seaSaveRastHard.RData")
#load("seaSaveRastEasy.RData")
#load("seaSaveEasomHard.RData")
#load("seaSaveEasomEasy.RData")
load("seaSaveEasomMed.RData")

isame = function(list1, list2){
        #Give theis function two vectors (ie. list1, list2). 
        #This function will determine if the vectors are the "same".
        #I define "same", as vectors containing the same number of elements in the same order.
        #This function returns a boolean if theye are or are not the same.
        truth = list1==list2
        n = length(truth)
        count=0
        for (item in truth){
                if (item){count=count+1}
        }

        return(n==count)
}


#typically keep will W*sLen (ie. number of points looking back (W) by the samples per point)
slideWindow = function(oldStuff, newStuff, keep){
        out = c(oldStuff, newStuff)
        outl = length(out)
        left = max((outl-(keep-1)), 1)
        out = out[left:outl]
        return( out )
}

test1 = function(int){
	#		 #use this for the two tailed test
	out = Z[int]>UCL #| Z[int]<LCL 
	out[1:ig1] = F
	return(out)
}



sig = 3#1#2#2.5
W = 10
#the larger this gets the more conservative the behaviour(less easy to flip bit) 
lambda = 0.35
flag = 1 
#
wSamples = c()
#
dev.new(width=14, height=14/3)
#
it = as.numeric(readline("Start Where? "))
ig1 = max(it-W, 1)
while(T){
	ZMax = seaSave[[1]][[it]]
	samples = seaSave[[2]][[it]]	
	maxes = seaSave[[3]][[it]]
	mcmcSize = seaSave[[4]][[it]]	
	#
	sEnd = length(samples)
        back = W*mcmcSize
        lEnd = max(1, sEnd-back)
        wSamples = samples[seq(lEnd, sEnd)]
	#
	wm = mean(wSamples)
        wv = var(wSamples)
        wmu = log( (wm^2)/(wv+wm^2)^0.5 )
        ws2 = log( 1 + (wv/(wm^2)) )
	#
	layout( matrix(c(1, 2, 3), nrow=1, ncol=3) )
	#Figure 1
	plot( seq(1, it), ZMax,
                "l",
                xlab="Sweeps",
                ylab="Objective Value",
                main="Best Z Value"
        )
	#
	it = it-1
	#
	UCL = matrix(NA, nrow=length(maxes), ncol=1)
	UCLN = matrix(NA, nrow=length(maxes), ncol=1)
	Z = matrix(NA, nrow=length(maxes), ncol=1)
	zMinus = 0	
	for( i in seq(1, it) ){
		flippi = it + 1 - i
		zMinus = lambda*maxes[flippi] + (1-lambda)*zMinus 
		Z[flippi] = zMinus
		#
		stuff = lambda/(2-lambda) * (1-(1-lambda)^(2*i))
        	UCL[flippi] = exp( wmu + (sig*(ws2*stuff)^0.5) )
		UCLN[flippi] = wm + (sig*(wv*stuff)^0.5)
	}
	#
	lmOut = btlm(X=seq(1, it), Z=Z[1:it], verb=0)	
	val = lmOut$trees[[length(lmOut$trees)]]$val
	valley = val[val!=0]
	print( val[val!=0] )
	if( flag ){
		ig1 = max(it-W, 1)
		if( (!isame(val, rep(0, length(val)))) ){
			ig1 = min( val[val>0] )
			flag = 0
		}
	} else if( length(valley)>1 ){
		ig1 = min(valley)
		W = it-max(valley)
	}	
	# 
	zMax = max(Z[1:it], UCL[1:it]+0.005, UCLN[1:it]+0.005)
	zMin = min(Z[1:it])
	#Figure 2
	plot( seq(1, it), Z[1:it],
		ylim=c(zMin, zMax),
		ylab="EWMA Value",
		xlab="Sweeps",
		main=sprintf("Sweeps: %d", it)
	)
	#	
	thresh = max(it-W, 1)
	out1 = test1(1:it)
	#
	points( seq(1, it)[out1], Z[1:it][out1],
                pch="O",
                col="red",
                cex=1.3
        )	
	#
	lines( seq(ig1, thresh),  UCL[ig1:thresh], col="red", lty=2 )	
	lines(  seq(thresh, it),   UCL[thresh:it], col="red" )	
	legend("topright",
		legend = c( sprintf("%1.1f Sigma log-N Interval", sig)),
		lty = c(1), 
		col = c("red")
	)
	#
	#Figure 3
	yMax = max(zMax, maxes[1:it])
        plot( seq(1, it), maxes[1:it],
                main="Maximum Expected Improvement",
                ylab='Max Improvement\n',
		ylim=c(zMin, yMax),
        )
	#segments(0, zMax+0.001, it+1, zMax+0.001, lty=3, cex=2)
	#
	move = readline("")
	if( move=="d" | move=="w" ){
		it = it+2
	}else if( move=="a" | move=="s" ){
		it = it
	}else{
		it = as.numeric(move) + 1
	}
	
}

