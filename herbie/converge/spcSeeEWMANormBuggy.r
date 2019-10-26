rm(list=ls())

library(qcc)

#load("seaSaveRoseEasyEasy.RData"); W=30
#load("seaSaveRoseEasy.RData"); W=10
#load("seaSaveRastHard.RData"); W=20
#load("seaSaveEasomHard.RData"); W=30
load("seaSaveShekelMed10.RData"); W=10

LNToN = function(LNSamples){
	wm = mean(LNSamples) #+ sqrt(10^-323)
	if(wm==0){
		print( sprintf("Check MCMC Convergence") )
		wm=NA
		wv=NA
		wmu=NA
		ws2=NA
		#           LN Stats    N Stats
		out = list( c(wm, wv), c(wmu, ws2) )
	}else{
        	wv = var(LNSamples)
        	#associated normal stats
        	wmu = log( (wm^2)/(wv+wm^2)^0.5 )
        	ws2 = log( 1 + (wv/(wm^2)) )
		#           LN Stats    N Stats
		out = list( c(wm, wv), c(wmu, ws2) )
	}
	
	return( out )
}

slice = function(x, n){
	out = split(x, as.integer((seq_along(x) - 1) / n))
	return( out )
}

getVect = function(outThing, LNNFlag, MVFlag){
	outVect = matrix(NA, nrow=length(outThing), ncol=1)
	if( LNNFlag=="LN" ){ lnn=1 }
	if( LNNFlag=="N" ){ lnn=2 }
	if( MVFlag=="M" ){ mv=1 }
	if( MVFlag=="V" ){ mv=2 }
	for( i in 1:length(outThing) ){
		outVect[i] = outThing[[i]][[lnn]][mv]
	}
	#
	return( outVect )
}

test1 = function(int){
        out = ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1] 
        #
        return( out )
}

wasd = function(){
	move = readline("")
	if( move=="d" | move=="w" ){
		it = it+1
	}else if( move=="a" | move=="s" ){
		it = it-1
	}else{
		it = as.numeric(move) + 1
	}
	return( it )
}


it0 = 30#W+1
it = it0
while(T){
	if( it!=it0 ){ dev.off() }
	#
	ZMax = seaSave[[1]][[it]]
	samples = seaSave[[2]][[it]]	
	#maxes = seaSave[[3]][[it]]
	mcmcSize = seaSave[[4]][[it]]	
	#
	sEnd = length(samples)
        back = W*mcmcSize
        lEnd = max(1, sEnd-back)
        wSamples = samples[seq(lEnd, sEnd)]
	mSamples = slice(samples, mcmcSize) 
	#
	wOut = LNToN(wSamples)
	wm = wOut[[1]][1]
	wv = wOut[[1]][2]
	wmu = wOut[[2]][1]
	ws2 = wOut[[2]][2]
	#
	mOut = lapply(mSamples, LNToN)
	mm = getVect(mOut, "LN", "M")
	mv = getVect(mOut, "LN", "V")
	mmu = getVect(mOut, "N", "M")
	ms2 = getVect(mOut, "N", "V")
	#
	brok = is.na(mmu)
	safe = !is.na(mmu)
	wSafe = brok[1:W]
	nsW = sum(wSafe)
	noW = sum(brok[W:length(safe)])
	#
	#Figure 1
	pdf("bestZ.pdf")
	plot( seq(1, it), ZMax,
                "l",
                xlab="Sweeps",
                ylab="Objective Value",
                main="Best Z Value"
        )
	dev.off()
	#
	#Figure 2 
	#think about how NA
	ewmaOut = ewma( rev(mmu[safe])[1:(W-nsW)], newdata=rev(mmu[safe])[(W-nsW):(length(mmu[safe])-noW)], plot=F )
	#
	testOut = rev(test1())
	dev.new()
	min2 = min(mmu[safe], ewmaOut$y, ewmaOut$limits)
	max2 = max(mmu[safe], ewmaOut$y, ewmaOut$limits)	
	plot( seq(1, it)[safe], rev(c(ewmaOut$statistics, ewmaOut$newstats)),
		pch=3, 
		ylim=c(min2, max2),
		xlab="Sweep Number",
		ylab="Summary Statistics",
		main="EWMA Convergence Chart"
	)
	print("hi")
	#lines( ewmaOut$x, rev(ewmaOut$y) )
	lines( seq(1, it)[safe], rev(ewmaOut$y) )
	lines( seq(1, it)[safe], rev(ewmaOut$limits[,1]), lty=2 )
	lines( seq(1, it)[safe], rev(ewmaOut$limits[,2]), lty=2 )
	#points( ewmaOut$x, rev(ewmaOut$y), pch=20 ) 
	points( seq(1, it)[safe], rev(ewmaOut$y), pch=20 )
	#points( ewmaOut$x[testOut], rev(ewmaOut$y)[testOut], pch=19, col="red" )		
	points( seq(1, it)[safe][testOut], rev(ewmaOut$y)[testOut], pch=19, col="red" )
	abline( v=(it+0.5)-W, lty=3 )
	abline( h=ewmaOut$center )

	#dev.copy2pdf(file="EWMA.pdf")		
	it = wasd()	
}








































#isame = function(list1, list2){
#        #Give theis function two vectors (ie. list1, list2). 
#        #This function will determine if the vectors are the "same".
#        #I define "same", as vectors containing the same number of elements in the same order.
#        #This function returns a boolean if theye are or are not the same.
#        truth = list1==list2
#        n = length(truth)
#        count=0
#        for (item in truth){
#                if (item){count=count+1}
#        }
#
#        return(n==count)
#}
#
##typically keep will W*sLen (ie. number of points looking back (W) by the samples per point)
#slideWindow = function(oldStuff, newStuff, keep){
#        out = c(oldStuff, newStuff)
#        outl = length(out)
#        left = max((outl-(keep-1)), 1)
#        out = out[left:outl]
#        return( out )
#}
#
#test1 = function(int){
#	#		 #use this for the two tailed test
#	out = Z[int]>UCL #| Z[int]<LCL 
#	#out[1:ig1] = F
#	return(out)
#}




#as.numeric(readline("Start Where? "))
#ig1 = 10
#sig = 2
#W = 30
#the larger this gets the more conservative the behaviour(less easy to flip bit) 
#lambda = 0.65
#flag = 1 
#
#wSamples = c()
#dev.new(width=12, height=14/3)
#dev.new(width=8, height=14/3)
#

#
	##Figure 1
	#layout( matrix(c(1, 2, 3, 4, 5, 6), nrow=2, ncol=3, byrow=T) )
        #plot( seq(1, it), ZMax,
        #        "l",
        #        xlab="Sweeps",
        #        ylab="Objective Value",
        #        main="Best Z Value"
        #)
        ##	
        #it = it-1
	#thresh = max(it-W, 1)
	##
	#LNSE = (mean(wv)/mcmcSize)^0.5
	#LNUCL = exp(wmu + sig*ws2^0.5)#LNSEws2^0.5
	#LNLCL = exp(wmu - sig*ws2^0.5)#LNSEws2^0.5
	##
	##Figure 2
	#min2 = min(mm, LNLCL)
	#max2 = max(mm, LNUCL)
	#plot( seq(1, it), mm,
	#	ylim=c(min2, max2)
	#)
	#segments(    ig1, LNUCL, thresh, LNUCL, col="red", lty=2)
	#segments( thresh, LNUCL,     it, LNUCL, col="red" )
	#segments(    ig1, LNLCL, thresh, LNLCL, col="red", lty=2)
	#segments( thresh, LNLCL,     it, LNLCL, col="red" )
	##
	#NSE = (mean(ws2)/mcmcSize)^0.5
	#NUCL = wmu + sig*ws2^0.5#NSE
	#NLCL = wmu - sig*ws2^0.5#NSE
	##
	##Figure 3
	#min3 = min(mmu, NLCL)
	#max3 = max(mmu, NUCL)
	#plot( seq(1, it), mmu,
	#	ylim=c(min3, max3)
	#)
	#segments(    ig1, NUCL, thresh, NUCL, col="red", lty=2)
	#segments( thresh, NUCL,     it, NUCL, col="red" )
	#segments(    ig1, NLCL, thresh, NLCL, col="red", lty=2)
	#segments( thresh, NLCL,     it, NLCL, col="red" )
	##
	##Figure 4
	#plot(1, 1)
	##
	#NEWMAUCL = matrix(NA, nrow=length(maxes), ncol=1)
	#NEWMALCL = matrix(NA, nrow=length(maxes), ncol=1)
	#NZ = matrix(NA, nrow=length(maxes), ncol=1)
	#LNEWMAUCL = matrix(NA, nrow=length(maxes), ncol=1)
	#LNEWMALCL = matrix(NA, nrow=length(maxes), ncol=1)
	#LNZ = matrix(NA, nrow=length(maxes), ncol=1)
	##
	#lnzMinus = 0
	#nzMinus = 0
	#for( i in 1:it ){
	#	flippi = it + 1 - i
	#	lnzMinus = lambda*mm[flippi] + (1-lambda)*lnzMinus
	#	nzMinus = lambda*mmu[flippi] + (1-lambda)*nzMinus
	#	LNZ[flippi] = lnzMinus
	#	NZ[flippi] = nzMinus
	#	#
	#	stuff = lambda/(2-lambda) * (1-(1-lambda)^(2*i))
	#	LNEWMAUCL[flippi] = exp( wmu + (sig*(ws2*stuff)^0.5) )
	#	LNEWMALCL[flippi] = exp( wmu - (sig*(ws2*stuff)^0.5) )
	#	NEWMAUCL[flippi] = wmu + (sig*((ws2)*stuff)^0.5)
	#	NEWMALCL[flippi] = wmu - (sig*((ws2)*stuff)^0.5)
	#}
	##
	##out1 = test1(1:it)
	##
	##Figure 5
	#min5 = min(LNZ[1:it], LNEWMALCL[1:it])
	#max5 = max(LNZ[1:it], LNEWMAUCL[1:it])
	#plot( 1:it, LNZ[1:it],
	#	ylim=c(min5, max5)
	#)
	#lines( ig1:thresh, LNEWMAUCL[ig1:thresh], col="red", lty=2 )
	#lines( ig1:thresh, LNEWMALCL[ig1:thresh], col="red", lty=2 )
	#lines(  thresh:it,  LNEWMAUCL[thresh:it], col="red", lty=1 )
	#lines(  thresh:it,  LNEWMALCL[thresh:it], col="red", lty=1 )
	##
	##Figure 6
        #min6 = min(NZ[1:it], NEWMALCL[1:it])
        #max6 = max(NZ[1:it], NEWMAUCL[1:it])
	#plot( 1:it, NZ[1:it],
	#	ylim=c(min6, max6)
	#)
	#lines( 1:thresh, NEWMAUCL[1:thresh], col="red", lty=2 )
	#lines( 1:thresh, NEWMALCL[1:thresh], col="red", lty=2 )
	#lines(  thresh:it,  NEWMAUCL[thresh:it], col="red", lty=1 )
	#lines(  thresh:it,  NEWMALCL[thresh:it], col="red", lty=1 )
	##

#
#	
#	#
#	wm = mean(wSamples)
#        wv = var(wSamples)
#        wmu = log( (wm^2)/(wv+wm^2)^0.5 )
#        ws2 = log( 1 + (wv/(wm^2)) )
#	#
#	layout( matrix(c(1, 2, 3), nrow=1, ncol=3) )
#	#Figure 1
#	plot( seq(1, it), ZMax,
#                "l",
#                xlab="Sweeps",
#                ylab="Objective Value",
#                main="Best Z Value"
#        )
#	#
#	it = it-1
#	#
#	UCL = matrix(NA, nrow=length(maxes), ncol=1)
#	LCL = matrix(NA, nrow=length(maxes), ncol=1)
#	Z = matrix(NA, nrow=length(maxes), ncol=1)
#	zMinus = 0	
#	for( i in seq(1, it) ){
#		flippi = it + 1 - i
#		zMinus = lambda*log(maxes[flippi]) + (1-lambda)*zMinus 
#		Z[flippi] = zMinus
#		#
#		stuff = lambda/(2-lambda) * (1-(1-lambda)^(2*i))
#        	UCL[flippi] = wm + (sig*(wv*stuff)^0.5)
#		LCL[flippi] = wm - (sig*(wv*stuff)^0.5)
#	}	
#	# 
#	zMax = max(Z[1:it], UCL[1:it]+0.01)
#	zMin = min(Z[1:it])
#	#Figure 2
#	plot( seq(1, it), Z[1:it],
#		#ylim=c(min(Z[1:it]), max(Z[1:it])),#c(zMin, zMax),
#		ylab="EWMA Value",
#		xlab="Sweeps",
#		main=sprintf("Sweeps: %d", it)
#	)
#	#	
#	thresh = max(it-W, 1)
#	#out1 = test1(1:it)
#	##
#	#points( seq(1, it)[out1], Z[1:it][out1],
#        #        pch="O",
#        #        col="red",
#        #        cex=1.3
#        #)	
#	##
#	#lines( seq(ig1, thresh),  UCL[ig1:thresh], col="red", lty=2 )	
#	#lines( seq(thresh, it),   UCL[thresh:it], col="red" )	
#	#lines( seq(ig1, thresh),  LCL[ig1:thresh], col="red", lty=2 )	
#	#lines( seq(thresh, it),   LCL[thresh:it], col="red" )
#	#legend("topright",
#	#	legend = c( sprintf("%1.1f Sigma log-N Interval", sig) ),
#	#	lty = c(1), 
#	#	col = c("red")
#	#)
#	#
#	wUpper = wmu + (sig*(ws2^0.5))
#        wLower = wmu - (sig*(ws2^0.5))
#	#
#	#Figure 3
#	min3 = min(wLower, lo)
#	yMax = max(zMax, maxes[1:it])
#        plot( seq(1, it), log(maxes[1:it]),
#                main="Maximum Expected Improvement",
#                ylab='Max Improvement\n',
#		ylim=c(min3, max3),
#        )
#	segments( thresh,  wUpper,     it, wUpper, col="red" )
#	segments( thresh,  wLower,     it, wLower, col="red" )
#	segments(    ig1,  wLower, thresh, wLower, col="red", lty=2 )
#	segments(    ig1,  wUpper, thresh, wUpper, col="red", lty=2 )

