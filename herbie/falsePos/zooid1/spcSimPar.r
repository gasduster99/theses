rm(list=ls())

#
library(tgp)
library(qcc)
#library(snow)
#library(Rmpi)
library(bigmemory)
suppressMessages(library(foreach, quietly=FALSE))
suppressMessages(library(doParallel, quietly=FALSE))

#
#FUNCTIONS
#

#
michal = function(xx, m=10, p=2){
  	##########################################################################
  	#
  	# INPUTS:
  	#
  	# xx = c(x1, x2)
  	# m = constant (optional), with default value 10
  	#
  	##########################################################################
  
	#
	xx = matrix(xx, ncol=p)
  	ii = 1:p
  	sum = rowSums(sin(xx) * (sin(sweep(xx^2/pi, 2, ii, '*')))^(2*m))
	#
	y = -sum
  	return(y)
}

#
grlee12 = function(xx){
	#
	xx = matrix(xx, ncol=1)
	#
	term1 = sin(10*pi*xx) / (2*xx)
	term2 = (xx-1)^4
  	#
	y = term1 + term2
  	return(y)
}

#
camel3 = function(xx){
	#
	xx = matrix(xx, ncol=2)
	x1 = xx[,1]
	x2 = xx[,2]
	#
	term1 = 2*x1^2
	term2 = -1.05*x1^4
	term3 = x1^6 / 6
	term4 = x1*x2
	term5 = x2^2
	#     
	y = term1 + term2 + term3 + term4 + term5
	return(y)
}

##
#levy = function(xx, p=2){
#	##########################################################################
#	#
#	# INPUT:
#	#
#	# xx = c(x1, x2, ..., xd)
#	#
#	##########################################################################
#	
#	#
#	xx = matrix(xx, ncol=p)
#	d = length(xx)
#	w = 1 + (xx - 1)/4
#	      
#	term1 = (sin(pi*w[1]))^2 
#	term3 = (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)
#	      
#	wi  = w[1:(d-1)]
#	sum = sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
#	      
#	y = term1 + sum + term3
#	return(y)
#}

#
mccorm = function(xx){
	#
	xx = matrix(x, ncol=2)
	x1 = xx[,1]
	x2 = xx[,2]
	#
	term1 = sin(x1 + x2)
	term2 = (x1 - x2)^2
	term3 = -1.5*x1
	term4 = 2.5*x2
	#     
	y = term1 + term2 + term3 + term4 + 1
	return(y)
}

#
rosenbrock = function(x){
        #
	x = matrix(x, ncol=2)
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

#
rastrigin = function(x, p=2){
        #
	x = matrix(x, ncol=p)
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
	if( wm==0 ){ wm=.Machine$double.eps }
	if( wv==0 ){ wv=.Machine$double.eps }
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
output = function(isOut, name, Zmax, it, mmu, lambda, ewma, W, test1, isSeeOut=T){
        #
        if(isOut){
                #Figure 1
                tMin = min(Zmax, na.rm=T)
                pdf(sprintf("bestZ%s.pdf", name))
                plot( seq(1, it), Zmax,
                        "l",
                        xlab="Iteration Number",
                        ylab="Objective Value",
                        main="Best Z Value",
                        ylim=c(tMin, max(Zmax, na.rm=T))
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
                min2 = min(mmu, ewmaOut$y, ewmaOut$limits, na.rm=T)
                max2 = max(mmu, ewmaOut$y, ewmaOut$limits, na.rm=T)
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
                if( isSeeOut ){ system(cmdString) }
		#
                system( sprintf('cp bestZ%s.pdf %s/pictures/', name, outPath) )
                system( sprintf('cp ewmaConvChart%s.pdf %s/pictures/', name, outPath) )
                system( sprintf('cp seeOut.pdf %s/pictures/seeOut%s.pdf', outPath, name) )
        }
}

#
monitorLog = function(isOut, rank, seed, it, itConv, grid){
	#look at output through:
	#~/.../falsePos$ while sleep <time>; do clear; echo -e "$(cat rank*/monitor.txt)"; done
	#or call bash monitor.sh which does that same thing 
	
	#
	every = 8
	if(isOut){
		#
		head = ''
		if(rank%%every==0){ head='[r]:[s] |' }	
		#
		line = sprintf('%03d:%03d |', rank, seed)
                #
		for(i in 1:length(grid)){
                        #
			if( rank%%every==0 ){ 
				#
				head = sprintf('%s w%03d', head, grid[i]) 
				#
				if( i==length(grid) ){ head=sprintf('%s\n', head) }
                        }
			#
			flag = '\\e[39m'
			if( unlist(itConv[i])==it ){ flag='\\e[92m' }
			line = sprintf('%s %s%4d', line, flag, unlist(itConv[i]))
                }	
		fBody = sprintf('%s%s\\e[39m', head, line)	

		#
		fc = file("monitor.txt")
		writeLines(fBody, fc)
		close(fc)
	}
}

#
meat = function(init, it){ #, dm){
        #
        #UNPACK INPUTS
        #

        #[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:mcmcSize; [[8]]:Xmax
        X = init[[1]]
        Z = init[[2]]
        Zmax = init[[3]]
        samples = init[[4]]
        maxes = init[[5]]
        out = init[[6]]
        sLen = init[[7]]
	Xmax = init[[8]]
        #nameCount = init[[8]]

        #
        #OPTIMIZE
        #

        #               
        out = suppressWarnings(tgp::optim.step.tgp( init$f, X=X, Z=Z, rect=rect, prev=out, improv=c(1,1), trace=T, verb=0, NN=NN ))
        ex = matrix(out$X, ncol=dm)
        fex = f(ex)
        X = rbind(X, ex)
        Z = c(Z, fex)
        Xmax = rbind(Xmax, X[Z==min(Z),])
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
                sLen=sLen,
		Xmax=Xmax
        ))
}

#
#MAIN
#

#fiddlers
M = 100 #8 #48 #100
threads = 46
NN = 200
makeOut = T
#
outPath = getwd() #'/home/nick/Documents/theses/herbie/falsePos/'
#
name = 'grlee12'
f = grlee12
rect = cbind(c(0.5), c(2.5))
opt = optim(0.5, grlee12, method='Brent', lower=0.5, upper=0.7)
zMin = opt$value
xMin = opt$par
wGrid = seq(20, 40, 2)
itMax = 200
##
#threads = 48
#name = 'mccormTry'
#f = mccorm
#rect = cbind(c(-1.5, -3), c(4, 4))
#zMin = -1.9133
#xMin = c(-0.54719, -1.54719)
#wGrid = seq(20, 40, 2)
#itMax = 200
##
#threads = 8
#name = 'rosenbrockNoCensor'
#f = rosenbrock
#rect = cbind(c(-2, -3), c(2, 5))
#zMin = 0
#xMin = c(1, 1)
#wGrid = seq(20, 40, 2)
#itMax = 200
##
#name = 'rastriginNoCensor'
#f = rastrigin
#rect = cbind(c(-2.5, -2.5), c(2.5, 2.5))
#zMin = 0
#xMin = rep(0, ncol(rect))
#wGrid = seq(20, 80, 5)
#itMax = 300
##

#
W = 40
lamGrid = seq(0.1, 0.65, 0.05)	#seq(0.1, 0.9, 0.05)
#xInitPerVol = 2
xInitPerD = 10
thresholdGrid = seq(4e-5, 1e-3, 4e-5) #seq(5e-5, 1e-3, 1e-4) #seq(2.75e-05, 7.75e-5, 5e-6)
threshold = 5e-4
#
rectVol = prod(rect[,2]-rect[,1])
dm = nrow(rect)

##
#not = c(17, 8, 20, 5, 23, 7)
#seed = 1:(M+length(not))
#seed = seed[!1:length(seed)%in%not]
seed = 1:M
#
pListDef = big.matrix(threads, 1,
	init           = -1, 
	backingfile    = 'pList.bin',
	descriptorfile = "pList.desc"
)
registerDoParallel(cores=threads)
out = foreach( m=1:M, .options.multicore=list(preschedule=F) )%dopar%{
#for( m in c(1) ){ #1:M ){
	#parallel
	pList = attach.big.matrix("pList.desc")
	pid = Sys.getpid()
	who = which(as.matrix(pList)==-1)
        if( m<=threads ){ pList[m] = pid
	} else{ pList[who[sample(1:length(who),1)]] = pid }
	rank = which(as.matrix(pList)==pid)-1 #(m-1)%%threads
	outName = sprintf('%s%s', name, m)
	#admin
	tic = Sys.time()	
	set.seed(seed[m])	
	#
	dir.create(sprintf('%s/rank%02d', outPath, rank))
	setwd(sprintf('%s/rank%02d', outPath, rank))
	system(sprintf('cp %s/seeOut.tex .', outPath))
	
	#sim init
	flags = as.data.frame( t(c(T, rep(F, length(thresholdGrid)), rep(F, length(wGrid)), rep(F, length(lamGrid)*length(wGrid)))) )
	nome = unlist(lapply(lamGrid, function(l) sprintf('ewmaL%.2fW%.2f', l, wGrid)))
	colnames(flags) = c('inLoop', sprintf('thresh%1.1e', thresholdGrid), sprintf('ewmaAutoW%.2f', wGrid), nome)	
	itConv = as.data.frame( t(c(rep(NA, length(thresholdGrid)), rep(NA, length(wGrid)), rep(NA, length(lamGrid)*length(wGrid)))) )
	colnames(itConv) = c(sprintf('thresh%1.1e', thresholdGrid), sprintf('ewmaAutoW%.2f', wGrid), nome)
	
	#optimization init
	xInit = xInitPerD*dm #Vol*rectVol
        X = lhs(xInit, rect)
        Z = f(X)
        Xmax = X[Z==min(Z),]
	Zmax = c(min(Z))
        samples = c()
        maxes = c()
        out = NULL
        sLen = NULL
	#[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:sLen; [[8]]:Xmax
	init = list(
		X=X, 
		Z=Z,
		Zmax=Zmax,
		samples=samples,
		maxes=maxes,
		out=out,
		sLen=sLen,
		Xmax=Xmax
	)
	#	
	autoLams = list()
	for(w in wGrid){ autoLams[[as.character(w)]]=numeric(0) }
	#
	it = 1
	while( flags$inLoop ){
		#
		init = meat(init, it)
		#[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:sLen; [[8]]:Xmax
		X = init[[1]]
		Z = init[[2]]
		Zmax = init[[3]]
		samples = init[[4]]
		maxes = init[[5]]
		out = init[[6]]
		sLen = init[[7]]
		Xmax = init[[8]]
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
		for(w in wGrid){
		
			#
			#Lambda Inference
			#
			
			#
			nome = sprintf('ewmaAutoW%.2f', w)
					
			#
			if( it>w ){
				#
                	        optOut = tryCatch( optimize(theFunk, c(0, 1), Zmax, samples, sLen, ewmaConvChart, w, slice, LNToN, getVect, ewma, ssError),
                	                                error=function(x){return(list(minimum=lambda))}
                	        )
                	        lambda = optOut$minimum
				autoLams[[as.character(w)]] = c(autoLams[[as.character(w)]], lambda)
				#
				ewmaOut = ewma( rev(mmu)[1:w], lambda=lambda, newdata=rev(mmu)[(w+1):length(mmu)], plot=F )
				#
                	        loggy = rev(ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1])
				place = length(loggy)-w
				loggyLeft = loggy[1:place]
                	        loggyRight = loggy[(place+1):length(loggy)]
				#
				flags[[nome]] = (any(loggyLeft) & !any(loggyRight)) | flags[[nome]]
				#if(!flags[[nome]]){ itConv[[nome]]=it+1 }
				#if(flags[[nome]] & it==(w+1)){ itConv[[nome]]=w+1 }
				#
				outNome = sprintf('%sEwmaLAutoW%.2f', outName, w)
				output( F, outNome, Zmax, it+1, mmu, lambda, ewma, w, test1 )
			}	
			
			#
			if(!flags[[nome]]){ itConv[[nome]]=it+1 }
			
			#
			#Lambda Grids
			#
			
			#
			for(lambda in lamGrid){
				#
				nome = sprintf('ewmaL%.2fW%.2f', lambda, w)
				#
				if( it>w ){
					#
					ewmaOut = ewma( rev(mmu)[1:w], lambda=lambda, newdata=rev(mmu)[(w+1):length(mmu)], plot=F )
					#
                		        loggy = rev(ewmaOut$y>ewmaOut$limits[,2] | ewmaOut$y<ewmaOut$limits[,1])
					place = length(loggy)-w
					loggyLeft = loggy[1:place]
                		        loggyRight = loggy[(place+1):length(loggy)]
                		        #
                		        flags[[nome]] = (any(loggyLeft) & !any(loggyRight)) | flags[[nome]]
                		        #if(!flags[[nome]]){ itConv[[nome]]=it+1 }
                		        #if(flags[[nome]] & it==(w+1)){ itConv[[nome]]=w+1 }	
					#
					outNome = sprintf('%sEwmaL%.2fW%.2f', outName, lambda, w)
					output( F, outNome, Zmax, it+1, mmu, lambda, ewma, w, test1, isSeeOut=F)
				}
				
				#
				if(!flags[[nome]]){ itConv[[nome]]=it+1 }
			}
		}
		
		#	
		monitorLog(makeOut, rank, seed[m], it+1, itConv[sprintf('ewmaAutoW%.2f', wGrid)], wGrid)
		
		#
		#Threshold
		#

		#
		for(threshold in thresholdGrid){
			nome = sprintf('thresh%1.1e', threshold)
                	flags[nome] = (min(mm)<threshold) | flags[nome]
                	if(!flags[nome]){ itConv[nome]=it+1 }
                	if(flags[nome] & it==1){ itConv[nome]=1 }
                }

		#
                it = it+1
                #if( (it-1)>W ){ output( makeOut, outName, Zmax, it, mmu, lambda, ewma, W, test1 ) }
		#itMax introduces right censoring.
                flags$inLoop = !all(unlist(flags[1+(1:length(thresholdGrid))])) & !all(unlist(flags[1+length(thresholdGrid)+(1:length(wGrid))]))   #!it>itMax #!( all(unlist(flags[-1])) | it>itMax )
	}
	#output structure
        outs = list(
                autoLams=autoLams,
                #X=X[(xInit+1):dim(X)[1],],
                #Xmax=Xmax[2:length(Zmax),],  
		#Zmax=Zmax[2:length(Zmax)],
                Xmax=Xmax[1:length(Zmax),],
		Zmax=Zmax[1:length(Zmax)], itConv=itConv, 
                origMean=mm,
                transMean=mmu,
                #loggyLeft=loggyLeft,
                #loggyRight=loggyRight,
		seed=seed[m],
		time=Sys.time()-tic
        )
	#
	pList[which(as.matrix(pList)==pid)] = -1
	#
	return(outs)
}

#
save.image(sprintf('%sL%.3f%.2fW%.2f%.2fM%d.RData', name, min(lamGrid), max(lamGrid), min(wGrid), max(wGrid), M))



#Grids:
#CHECK lambda; autoLam
#Ws
#thesholds
#xInitPerVol
#seed
#run time


#
#JUNK
#

#flags = data.frame(
	#	inLoop = T,
	#	thresh = F, 
	#	ewma   = F	
	#)
#itConv = data.frame(
	#	thresh = NA,
	#	ewma   = NA
	#)
