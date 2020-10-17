rm(list=ls())

#
library(tgp)
library(qcc)
#library(snow)
#library(Rmpi)
suppressMessages(library(foreach, quietly=FALSE))
suppressMessages(library(doParallel, quietly=FALSE))

#
#FUNCTIONS
#

#
rosenbrock = function(x){
        #
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

#
rastrigin = function(x){
        #
        out = matrix(NaN, nrow=dim(x)[1], ncol=1)
        for (i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = 10*length(ex)+sum(ex^2-10*cos(2*pi*ex))
        }
        return(out)
}

#
easom = function(x){
        #out = -apply(cos(x), 1, prod) * exp(-apply((x-pi)^2, 1, sum))
        out = matrix(NaN, nrow=dim(x)[1], ncol=1)
        for(i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = -prod(cos(ex)) * exp( -sum((ex-pi)^2) )
        }
        return(out)
}

#
test1 = function(ewmaOut){
        #
        out = ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1]
        #
        return(out)
}

#typically keep will be W*sLen (ie. number of points looking back (W) by the samples per point)
slideWindow = function(oldStuff, newStuff, keep){
        #
        out = c(oldStuff, newStuff)
        outl = length(out)
        left = max((outl-(keep-1)), 1)
        out = out[left:outl]
        #
        return( out )
}

#
LNToN = function(LNSamples){
        #log-normal stats
        wm = mean(LNSamples)
        wv = var(LNSamples)
        #associated normal stats
        wmu = log( (wm^2)/(wv+wm^2)^0.5 )
        ws2 = log( 1 + (wv/(wm^2)) )
        #
        #           LN Stats    N Stats
        out = list( c(wm, wv), c(wmu, ws2) )
}

#
slice = function(x, n){
        #
        out = split(x, as.integer((seq_along(x) - 1) / n))
        return(out)
}

#
getVect = function(outThing, LNNFlag, MVFlag){
        #
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

#
ewmaConvChart = function(lamb, Zmax, samples, mcmcSize, W, slice, LNToN, getVect, ewma){
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
        ewmaOut = ewma( rev(mmu)[1:W], lambda=lamb, newdata=rev(mmu)[(W):length(mmu)], plot=F )
        testOut = rev( ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1] )
        #fudge = sample(rev(sort(ewmaOut$statistics))[1:3], 1)
        min2 = min(mmu, ewmaOut$y, ewmaOut$limits)
        max2 = max(mmu, ewmaOut$y, ewmaOut$limits)
        #
        out = list(
                rev(c(ewmaOut$statistics, ewmaOut$newstats)),   #1 z
                ewmaOut$x[!testOut],                            #2
                rev(ewmaOut$y)[!testOut],                       #3
                ewmaOut$x[testOut],                             #4
                rev(ewmaOut$y)[testOut],                        #5
                ewmaOut$x,                                      #6
                rev(ewmaOut$y),                                 #7 \hat{z}
                matrix(rev(ewmaOut$limits), ncol=2),            #8
                lamb                                            #9
        )
        return( out )
}

#
ssError = function(ewmaOut){
        #
        lambda = ewmaOut[[9]]
        en = length(ewmaOut[[7]])
        zHat = ewmaOut[[1]][1:(en-1)]
        z = ewmaOut[[7]]
        #
        z1Hat = z[1:en-1]*lambda + zHat*(1-lambda)
        #
        return( sum((z1Hat-z[2:en])^2) )
}

#
theFunk = function(lamb, Zmax, samples, mcmcSize, ewmaConvChart, W, slice, LNToN, getVect, ewma, ssError){
        #
        ewmaOut = ewmaConvChart(lamb, Zmax, samples, mcmcSize, W, slice, LNToN, getVect, ewma)
        return( ssError(ewmaOut) )
}

#
output = function(isOut, name, Zmax, it, mmu, lambda, ewma, W, test1){
        #
        if(isOut){
                #Figure 1
                tMin = min(Zmax)
                pdf(sprintf("bestZ%s.pdf", name))
                plot( seq(1, it), Zmax,
                        "l",
                        xlab="Iteration Number",
                        ylab="Objective Value",
                        main="Best Z Value",
                        ylim=c(tMin, max(Zmax))
                )
                #abline( h=tMin, lty=3 )
                #legend("topright", c("\nBest\n Z Value", "\nLiterature \nMinimum", NA),
                #       lty=c(1, 3, NA)
                #)
                dev.off()
                #
                ewmaOut = ewma( rev(mmu)[1:W], lambda=lambda, newdata=rev(mmu)[(W):length(mmu)], plot=F )
                testOut = rev(test1(ewmaOut))
                #fudge = sample(rev(sort(ewmaOut$statistics))[1:3], 1)#sort(mm)[1:5], 1)
                min2 = min(mmu, ewmaOut$y, ewmaOut$limits)
                max2 = max(mmu, ewmaOut$y, ewmaOut$limits)
                #
                #Figure 2
                #dev.new()
                pdf(sprintf("ewmaConvChart%s.pdf", name))
                plot( seq(1, it), rev(c(ewmaOut$statistics, ewmaOut$newstats)),
                        pch=3,
                        ylim=c(min2, max2+1),
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
                cmdString = sprintf("pdflatex '\\def\\bestName{bestZ%s.pdf} \\def\\ewmaName{ewmaConvChart%s.pdf} \\input{seeOut.tex}'", name, name)
                system( cmdString )
                #system( sprintf('mv bestZ%s.pdf /home/nick/Documents/falsePos/pictures/', name) )
                #system( sprintf('mv ewmaConvChart%s.pdf /home/nick/Documents/falsePos/pictures/', name) )
                #system( sprintf('cp seeOut.pdf /home/nick/Documents/falsePos/pictures/seeOut%s.pdf', name) )
                system( sprintf('mv bestZ%s.pdf %s/pictures/', name, outPath) )
                system( sprintf('mv ewmaConvChart%s.pdf %s/pictures/', name, outPath) )
                system( sprintf('cp seeOut.pdf %s/pictures/seeOut%s.pdf', outPath, name) )
        }
}

#
meat = function(init, it){ #, dm){
        #
        #UNPACK INPUTS
        #

        #[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:mcmcSize
        X = init[[1]]
        Z = init[[2]]
        Zmax = init[[3]]
        samples = init[[4]]
        maxes = init[[5]]
        out = init[[6]]
        sLen = init[[7]]
        #nameCount = init[[8]]

        #
        #OPTIMIZE
        #

        #               
        out = tgp::optim.step.tgp( init$f, X=X, Z=Z, rect=rect, prev=out, improv=c(1,1), trace=T, verb=0, NN=NN )
        ex = matrix(out$X, ncol=dm)
        fex = f(ex)
        X = rbind(X, ex)
        Z = c(Z, fex)
        Zmax = c(Zmax, min(Z))
        #
        improvSamples = out$obj$trace$preds$improv
        EimprovAll = out$obj$improv
        maxI = which( EimprovAll$rank==1 )
        maxSamples = unlist( improvSamples[maxI] )
        samples = c(samples, maxSamples)
        #
        m = mean(maxSamples)
        maxes = c(maxes, m)
        sLen = length(maxSamples)

        #
        return(list(
		X=X,
                Z=Z,
                Zmax=Zmax,
                samples=samples,
                maxes=maxes,
                out=out,
                sLen=sLen
        ))
}

#
#MAIN
#

#fiddlers
M = 8 #1000
threads = 8
NN = 200
makeOut = T
#
outPath = getwd() #'/home/nick/Documents/theses/herbie/falsePos/'
name = 'test'
f = rosenbrock
rect = cbind(c(-2, -3), c(2, 5))
dm = nrow(rect)
W = 30

#
threshold = 5e-4

#
registerDoParallel(cores=threads)
out = foreach( m=1:M )%dopar%{
#for( m in 1:M ){
	#parallel
	rank = (m-1)%%threads
	outName = sprintf('%s%s', name, m)
	#
	dir.create(sprintf('%s/rank%s', outPath, rank))
	setwd(sprintf('%s/rank%s', outPath, rank))
	system(sprintf('cp %s/seeOut.tex .', outPath))
	
	#sim init
	flags = data.frame(
		inLoop = T,
		thresh = F, 
		ewma   = F	
	)		
	itConv = data.frame(
		thresh = NA,
		ewma   = NA
	)
	it = 1
	#optimization init
	#dm = nrow(rect) 
        X = lhs(40, rect)
        Z = f(X)
        Zmax = c(min(Z))
        samples = c()
        maxes = c()
        out = NULL
        sLen = NULL
	#[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:sLen; 
	init = list(
		X=X, 
		Z=Z,
		Zmax=Zmax,
		samples=samples,
		maxes=maxes,
		out=out,
		sLen=sLen
	)
	#
	while( flags$inLoop ){
		#
		init = meat(init, it)
		#[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:sLen; 
		X = init[[1]]
		Z = init[[2]]
		Zmax = init[[3]]
		samples = init[[4]]
		maxes = init[[5]]
		out = init[[6]]
		sLen = init[[7]]	
		#
		sEnd = length(samples)
		back = W*sLen
		lEnd = max(1, sEnd-back)
		wSamples = samples[seq(lEnd, sEnd)]
		mSamples = slice(samples, sLen)
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
		if( it>W ){
			 #
                        optOut = tryCatch( optimize(theFunk, c(0, 1), Zmax, samples, sLen, ewmaConvChart, W, slice, LNToN, getVect, ewma, ssError),
                                                error=function(x){return(list(minimum=lambda))}
                        )
                        lambda = optOut$minimum
                        #optOut = optim(0.4, theFunk, Zmax, samples, sLen, ewmaConvChart, W, slice, LNToN, getVect, ewma, ssError)
                        #lambda = optOut$par
                        ewmaOut = ewma( rev(mmu)[1:W], lambda=lambda, newdata=rev(mmu)[(W+1):length(mmu)], plot=F )
                        #
                        loggy = rev(ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1])
                        place = length(loggy)-W
                        loggyLeft = loggy[1:place]
                        loggyRight = loggy[(place+1):length(loggy)]
                        flags$ewma = (any(loggyLeft) & !any(loggyRight)) | flags$ewma
                        #
                        flags$inLoop = flags$ewma & flags$thresh
                        #
                        lambdas = c(lambdas, lambda)
                        if(!flags$ewma){ itConv$ewma=it+1 }
                        if(flags$ewma & it==(W+1)){ itConv$ewma=W+1 }	
		}
		#
                flags$thresh = (min(mm)<threshold) | flags$thresh
                if(!flags$thresh){ itConv$thresh=it+1 }
                if(flags$thresh & it==1){ itConv$thresh=1 }
                #
                it = it+1
                if( (it-1)>W ){ output( makeOut, nameCount, Zmax, it, mmu, lambda, ewma, W, test1 ) }
	}
	#output structure
        outs = list(
                lambdas=lambdas,
                X=X[41:dim(X)[1],],
                Zmax=Zmax[2:length(Zmax)],
                itConv=itConv, 
                origMean=mm,
                transMean=mmu,
                loggyLeft=loggyLeft,
                loggyRight=loggyRight
        )
	#print(m)    
	return(outs)
}


#grid of lambdas
#grid of Ws
#
