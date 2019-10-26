rm(list=ls())

library(qcc)

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

wasd = function(name){
	cmdString = sprintf("pdflatex '\\def\\bestName{ss%s.pdf} \\def\\ewmaName{ewmaConvChart%s.pdf} \\input{seeOut.tex}'", name, name)
	writeLines( cmdString )
	system( cmdString )
	move = readline("")
	if( move=="d" | move=="w" ){
		it = it+1
	}else if( move=="a" | move=="s" ){
		it = it-1
	}else if( move=="save" | move=="S" | move=="Save" | move=="SAVE" ){
		add = readline("Distinguishing Text? ")
		#saveBest = sprintf("cp ./bestZ%s.pdf ../../../../../../nickHerbie/figures/bestZ%s%s.pdf", name, name, add)
		saveSS = sprintf("cp ./ss%s.pdf ../../../../../../nickHerbie/figures/ss%s%s.pdf", name, name, add)
		saveEWMA = sprintf("cp ./ewmaConvChart%s.pdf ../../../../../../nickHerbie/figures/ewmaConvChart%s%s.pdf", name, name, add)
		system( saveSS )
		system( saveEWMA )
	}else{
		it = as.numeric(move) + 1
	}	
	return( it )
}

ewmaConvChart = function(lamb, it){
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
	ewmaOut = ewma( rev(mmu)[1:W], lambda=lamb, newdata=rev(mmu)[W:length(mmu)], plot=F )
        testOut = rev( ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1] )
        fudge = sample(rev(sort(ewmaOut$statistics))[1:3], 1)
        min2 = min(mmu, ewmaOut$y, ewmaOut$limits)
        max2 = max(mmu, ewmaOut$y, ewmaOut$limits)
	#
	out = list(
		seq(1, it), 					#1
		rev(c(ewmaOut$statistics, ewmaOut$newstats)), 	#2 z
		ewmaOut$x[!testOut],				#3
		rev(ewmaOut$y)[!testOut],			#4
		ewmaOut$x[testOut],				#5
		rev(ewmaOut$y)[testOut],			#6
		ewmaOut$x,					#7
		rev(ewmaOut$y),					#8 \hat{z}
		matrix(rev(ewmaOut$limits), ncol=2),		#9
		lamb						#10
	)
	return( out )
}

ewmaPlot = function(ewmaOut, it){
	plot(ewmaOut[[1]], ewmaOut[[2]], 
		pch=3,
		#ylim=c(min2, max2),
               	xlab="Iteration Number",
               	ylab="Summary Statistics",
               	main="EWMA Convergence Chart"
	)
	points(ewmaOut[[3]], ewmaOut[[4]], pch=20)
	points(ewmaOut[[5]], ewmaOut[[6]], col='red')
	lines(ewmaOut[[7]], ewmaOut[[8]])
	lines(ewmaOut[[1]], ewmaOut[[9]][,1], lty=2)
	lines(ewmaOut[[1]], ewmaOut[[9]][,2], lty=2)
	abline(v=(it+0.5)-W, lty=3)
}

ssError = function(ewmaOut){
	lambda = ewmaOut[[10]] 
	en = length(ewmaOut[[8]])
	zHat = ewmaOut[[2]][1:(en-1)]
	z = ewmaOut[[8]] 
	#
	z1Hat = z[1:en-1]*lambda + zHat*(1-lambda)
	#
	return( sum((z1Hat-z[2:en])^2) )
}

findEWMA = function(it, h){
	ls = seq(0.01, 1, h)
	nl = length(ls)
	ils = seq(1, nl)
	ss = matrix(NA, nrow=nl, ncol=1)
	ewma = NA
	for(i in ils){
		ewmaLast = ewma
	        ewma = ewmaConvChart(ls[i], it)
	        ss[i] = ssError(ewma)
	        #if( ss[i]<=min(ss[1:i]) ){ ewmaOut=ewma }
		if( ss[i]>min(ss[1:i]) ){ return(ewmaLast) }
	}
	#return( ewmaOut )
}

findSS = function(it, h){
	ls = seq(0.01, 1, h)
	nl = length(ls)
	ils = seq(1, nl)
	ss = matrix(NA, nrow=nl, ncol=1)
	for(i in ils){
	        ewma = ewmaConvChart(ls[i], it)
	        ss[i] = ssError(ewma)
	        #if( ss[i]<=min(ss[1:i]) ){ ewmaOut=ewma }
	}
	return( list(ls, ss) )
}

#it0 = W+1
#it = it0
#while(T){
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































