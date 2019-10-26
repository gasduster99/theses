rm(list=ls())

library(qcc)

#name = "RosegpPic"; uu=30; tMin=0; lambda=0.2
name = "RoseEasyEasy"; uu=30; tMin=0; lambda=0.2
#
#name = "RastHard"; uu=60; tMin=0; lambda=0.3
##name = "RastHard2"; uu=30; tMin=0; lambda=0.3

#name = "EasomMed"; uu=20; tMin=-1; lambda=0.1


##name = "RastHard"; uu=20; tMin=0; lambda=0.2
#name = "EasomMedHard"; uu=15; tMin=-1; lambda=0.02
#name = "EasomHard"; uu=30; tMin=-1; lambda=0.3
#load("seaSaveRoseEasyEasy.RData"); W=30
#load("seaSaveRoseEasy.RData"); W=10
#load("seaSaveRastHard.RData"); W=20
#load("seaSaveEasomHard.RData"); W=30
#load("seaSaveShekelMed10.RData"); W=20

load(sprintf("seaSave%s.RData", name))
W = uu

fudge = 0
LNToN = function(LNSamples){
	wm = mean(LNSamples)
	#a little fudggy goodness 
	if (wm==0){ wm=wm+fudge }
        wv = var(LNSamples)
        #associated normal stats
        wmu = log( (wm^2)/(wv+wm^2)^0.5 )
        ws2 = log( 1 + (wv/(wm^2)) )
	#
	#           LN Stats    N Stats
	out = list( c(wm, wv), c(wmu, ws2) )
}

slice = function(x, n){
	out = split(x, as.integer((seq_along(x) - 1) / n))
	return(out)
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
        return(out)
}

wasd = function(name){
	cmdString = sprintf("pdflatex '\\def\\bestName{bestZ%s.pdf} \\def\\ewmaName{ewmaConvChart%s.pdf} \\input{seeOut.tex}'", name, name)
	writeLines( cmdString )
	system( cmdString )
	move = readline("")
	if( move=="d" | move=="w" ){
		it = it+1
	}else if( move=="a" | move=="s" ){
		it = it-1
	}else if( move=="save" | move=="S" | move=="Save" | move=="SAVE" ){
		add = readline("Distinguishing Text? ")
		saveBest = sprintf("cp ./bestZ%s.pdf ../../nickHerbie/figures/bestZ%s%s.pdf", name, name, add)
		saveEWMA = sprintf("cp ./ewmaConvChart%s.pdf ../../figures/ewmaConvChart%s%s.pdf", name, name, add)
		system( saveBest )
		system( saveEWMA )
	}else{
		it = as.numeric(move) + 1
	}	
	return(it)
}


it0 = W+1
it = it0
while(T){
	#if( it!=it0 ){ dev.off() }
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
	#Figure 1
	pdf(sprintf("bestZ%s.pdf", name))
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
	ewmaOut = ewma( rev(mmu)[1:W], lambda=lambda, newdata=rev(mmu)[W:length(mmu)], plot=F )	
	testOut = rev(test1())
	fudge = sample(rev(sort(ewmaOut$statistics))[1:3], 1)
	min2 = min(mmu, ewmaOut$y, ewmaOut$limits)
	max2 = max(mmu, ewmaOut$y, ewmaOut$limits)
	#
	#Figure 2
	#dev.new()
	pdf(sprintf("ewmaConvChart%s.pdf", name))
	plot( seq(1, it), rev(c(ewmaOut$statistics, ewmaOut$newstats)),
		pch=3, 
		ylim=c(min2, max2),
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
	legend( "topright", c(expression(Y[i]), expression(Z[i]*" Violation"), expression(Z[i]*" Control"), "Control Limit", NA,"\nControl Window\nBoundary"), 
		pch=c(3 , 19, 20, NA, NA, NA),
		lty=c(NA, NA, NA, 2, 3, NA),
		col=c("black", "red", "black", "black", "black", NA)
	)
	dev.off()
	
	#
	#it = it-1
	##	
	#if( it>it0 ){ W=ceiling(W/V*var(c(ewmaOut$statistics, ewmaOut$newstats))) }
	#V = var(c(ewmaOut$statistics, ewmaOut$newstats))
	#print(W)
	##
	it = wasd(name)	
}































