rm(list=ls())

#
library(VGAM)
library(boot)
library(Ryacas)
library(pracma)
library(mvtnorm)
library(graphics)
library(parallel)
library(latex2exp)
library(rootSolve)
library(matrixStats)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)
#
source('ddClass0.1.1.r')
#source('gpFunk.r')
source('gpFunkGaussLxLy.r')

#
#FUNCTIONS
#

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
SRR = function(B, alpha, beta, gamma){
        alpha*B*(1-beta*gamma*B)^(1/gamma)
}

#Equilibrium Equations

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
Bbar_r = as_r( with_value(Bbar, "F", ysym("FF")) )
BBar = function(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#
FBbar = ysym("F*B")
FBbar = with_value(FBbar, "B", Bbar)
FBbar = with_value(FBbar, "F", ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
dFBdF_r = as_r(dFBdF)
FDebug = function(FF, M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        eval(dFBdFexpr)
}
FMsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        #
	out <- tryCatch({
		uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10))$root
	}, error=function(err){
                out = NA
        })
	#
	return(out)
}

#
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#RP Inversion Tools

#beta does not matter for either of getAlphaFmsy or getGammaFmsy
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma, dFBdFexpr=dFBdF_r){
        #
        #capture.output(
        out <- tryCatch({
                uniroot(function(alpha){ eval(dFBdFexpr) }, c(eps(), 100))$root
        }, error=function(err){
                out = NA
        })#, file="/dev/null")
        #
        return(out)
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta, dFBdFexpr=dFBdF_r){
        #FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(dFBdFexpr) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
        #
	f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
	#
	out <- tryCatch({
        	uniroot(f, c(0, 10), tol=eps())$root
	}, error=function(err){
		out = NA
	})
	return(out)
}

#gamma|alpha, zeta
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
        #
	f = function(g){
                BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - zeta
        }
	out <- tryCatch({
        	uniroot(f, c(-10, 10))$root
	}, error=function(err){
                out = NA
        })
        return(out)
}

#
getZeta = function(FF, M, k, w, W, alpha, beta, gamma){
        BBar(FF, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
}

#
getZetaBH = function(x, M, k, W, aS, a0){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
        #
        gamma = -1
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaBH = Vectorize(getZetaBH, "x")

##
#getZetaBH = function(x, M, k, W, a0){
#	#
#        if( is.na(x) ){ return(NA) }
#	#
#        w = W*(1-exp(-k*a0))
#        #
#        gamma = -1
#        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
#        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
#        #
#        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
#        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
#        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BZero #BBar(0, M, k, w, W, alpha, beta, gamma)
#        #
#        return(zetaHat)
#}
#getZetaBH = Vectorize(getZetaBH, "x")

#
invert = function(zeta, xi, B0, M, k, w, W, alphaStart=1, betaStart=1, gammaStart=1){
        #
        alpha = alphaStart
        beta  = betaStart
        gamma = gammaStart
        #
        xTol = 10^-3
        zTol = xTol
        bTol = 10
        xiHat = xi+xTol*10
        zetaHat = zeta+zTol*10
        #
        while( abs(xiHat-xi)>=xTol | abs(zetaHat-zeta)>=zTol ){
                #
                gamma = getGammaZeta(zeta, xi*M, M, k, w, W, alpha, beta)
                alpha = getAlphaFmsy(xi*M, M, k, w, W, beta, gamma)
                beta  = getBeta(B0, M, k, w, W, alpha, gamma)
                #
                BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
                xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
                zetaHat = BBar(xiHat*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
                #
                #xis = c(xis, xiHat)
                #zetas = c(zetas, zetaHat)
                #i = i+1
        }
        #
        return( data.frame(alpha=alpha, beta=beta, gamma=gamma) )
}

#Statistics

#
myDist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-getZetaBH(x, M, kappa, WW, aS, a0))^2) }

#
getlFV = function(fit, MM=10^4, samples=F){
	#
	ww = vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0) #fit$WW*(1-exp(-fit$kappa*fit$a0))
	#
	who = c('lalpha')
	C = fit$rsCov[who, who]
	m = c(fit$lalpha)
	sam = rnorm(MM, m, sqrt(C)) #rmvnorm(M, m, C)
	#FMsy(a, fit$gamma, M))})
	als = exp(sam)
	lFs = sapply(als, function(a){log(FMsy(fit$M, fit$kappa, ww, fit$WW, a, fit$beta, fit$gamma))}) 
	#
	out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T))
	if( samples ){ out$lFSamples=lFs }
	#
	return(out)
}

#
getlV = function(fit, MM=10^4, samples=F){
	#
	ww = vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0) #fit$WW*(1-exp(-fit$kappa*fit$a0))
	#
	who = c('lalpha', 'lbeta')
	C = fit$rsCov[who, who]
	m = c(fit$lalpha, fit$lbeta)
	sam = rmvnorm(MM, m, C) #rnorm(MM, m, sqrt(C)) #
	#NOTE: modify FMsy and BBar
	abs = exp(sam)
	lFs = sapply(abs[,1], function(a){log(FMsy(fit$M, fit$kappa, ww, fit$WW, a, fit$beta, fit$gamma))})
	#lKs = mapply(function(p1, p2){log(PBar(p1, p2, fit$gamma, 0, M))}, abs[,1], abs[,2])
	#(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r)
	lB0s = mapply(function(p1, p2){log(BBar(0, fit$M, fit$kappa, ww, fit$WW, p1, p2, fit$gamma))}, abs[,1], abs[,2])
	#
	out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T), lB0=mean(lB0s, na.rm=T), lB0V=var(lB0s, na.rm=T))
	if( samples ){ out$lFSamples=lFs; out$lB0Samples=lB0s }
	#
	return(out)
}

#
getData = function(dir, xiRange, zetaRange){
	#
	
	#
	fitFiles = sprintf( "%s%s", dir, list.files(path=dir, pattern=glob2rx("fit_*.rda")) )
	#
	i = 1
	D = data.frame(
		xiBin=double(), 
		zetaBin=double(), 
		xiSeed=double(), 
		zetaSeed=double(), 
		xiHat=double(), 
		zetaHat=double(),
		lF=double(), 
		lFV=double(), 
		lB0=double(), 
		lB0V=double(),
		minDist=double(), 
		stringsAsFactors=F
	)
	#
	for(f in fitFiles){
		#
		#print(f)
		fit = readRDS(f)
		dat = readRDS(gsub("fit", "datGen", f))			

		#	
		if( dat$xi<xiRange[1] | dat$xi>xiRange[2] ){
			#| dat$zeta<zetaRange[1] | dat$zeta>zetaRange[2]){ 
			next 
		}
		
		#
		xiBin   = strsplit(f, "_")[[1]]
		zetaBin = xiBin[3]
		xiBin   = xiBin[2]
		#
		zetaBin = gsub("zeta", "", zetaBin)
		zetaBin = as.numeric(gsub(".rda", "", zetaBin))	
		xiBin = as.numeric(gsub("xi", "", xiBin))

		#
		#Fs = FMsy(fit$alpha, fit$gamma, M)
		ww = vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0) #fit$WW*(1-exp(-fit$kappa*fit$a0))
		Fs = FMsy(fit$M, fit$kappa, ww, fit$WW, fit$alpha, fit$beta, fit$gamma)
		xiHat = Fs/fit$M
		#
		#BStar = PBar(fit$alpha, fit$beta, fit$gamma, Fs, M)		
		BStar = BBar(Fs, fit$M, fit$kappa, ww, fit$WW, fit$alpha, fit$beta, fit$gamma)
		#BZero = PBar(fit$alpha, fit$beta, fit$gamma, 0, M)
		BZero = BBar(0, fit$M, fit$kappa, ww, fit$WW, fit$alpha, fit$beta, fit$gamma)
		zetaHat = BStar/BZero	
		
		#
		md = stats::optimize(myDist, c(0, max(xiRange)), xi=dat$xi, zeta=dat$zeta)$objective
		
		#
		if(length(fit$rsCov)==0){ lfv=0; lbv=0 }else{ 
			#
			lV = getlV(fit)
			#
			lfv = lV$lFV
			lbv = lV$lB0V 
		}

		#
		D[i,] = c(xiBin, zetaBin, dat$xi, dat$zeta, xiHat, zetaHat, log(Fs), lfv, log(BZero), lbv, md)	
		i = i+1
	}
	#
	return(D)
}

#
edge = function(HS){
	#
	line = HS
	for(row in 1:nrow(HS)){
	        l = HS[row,]
	        who = which(l)[1]
	        line[row,] = rep(F, length(l))
	        line[row,who] = T
	}
	#
	return(line)
}

#
edgeBH = function(HS){
	#
	line = HS
	for(row in 1:nrow(HS)){
	        #
		zetaBH = bh(xiStar[row])
		#zetaOne = which(zetaStar>zetaBH)[1]
		#
		l = HS[row,] & zetaStar>zetaBH
	        who = which(l)[1]
	        line[row,] = rep(F, length(l))
	        line[row,who] = T
	}
	#
	return(line)
}

#
hsToLine = function(x, HS){
	#
	lineMat = edge(HS)
	bound = c(bound, rev( ksmooth(x, ksmooth(x,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y ))

}


##
##DATA
##
#
##
##mod = #"ExpT45N150K1" #"FlatT30N150A15K0.1" #"ExpT45N150A15" #"ExpT45N150A7.5" # #"ExpT45N150Wide"
##fail
##mod = "FlatT30N150A-1AS15K0.1"
##fail
##mod = "FlatT30N150A-0.5AS15K0.1" 
###worked w/ F help w/o K help
##mod = "ExpT45N150A-1AS2" 
###worked w/ F help w/o K help
##mod = "ExpT45N150A-1AS15K0.1" 
##fail
##mod = "ExpT45N150A-0.5AS15K0.1" 
#
### worked but no xi<1
##mod = "ExpT45N300AS1K1" 
##kinda worked
##mod = "ExpT45N300AS10K10"
#
###failed
##mod = "ExpT45N300AS10K0.1"
##WORKED
##mod = "ExpT45N300A0-1AS10K0.1"; fv=1
##GOOD START: add a bit of refinement
##mod = "ExpT45N150A0-1AS2K0.1"
##failed
##mod = "ExpT45N300A0-1AS10K0.1N28"
###WORKED 
##mod = "ExpT45N300AS0.1K10"; fv=1
##
##mod = "ExpT45N300A0-1AS0.1K0.1"; fv=10
##
##mod = "ExpT45N150A0-1AS4K0.2"; fv=1
##mod = "ExpT45N150A0-1AS4K0.2N28"; fv=1
#
#
###failed Fmsy overestimated
##mod = "FlatT45N300A0-1AS10K0.1"
##bumpy 
##mod = "FlatT45N300A0-1AS10K0.1N56"
##mod = "FlatT45N300A0-1AS10K0.1N84"
#
###GOOD START: add refinement
##mod = "FlatT45N150A0-1AS0.1K10N56"; fv=100
#
###GOOD START: add refinment
##mod = "FlatT45N150A0-1AS2K0.1N84Edge"; fv=200
#
###
##mod = "FlatT45N150A0-1AS4K0.2N56"; fv=100
#
##
#mod = "FlatT45N150A0-1AS1K0.5N56"; fv=100
#
##
#place = sprintf("./modsDD%s/", mod)
#
##
#xiRes = 0.5
#zetaTop = 0.6 #0.7
#zetaBot = 0.2 #0.1
#xiBot = 0.5
#xiTop = 3.5
#
##
#i = 1
#C = NULL
#while( is.null(C) ){
#	#
#	f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit_*.rda"))[i] )
#	fOne = readRDS(f)
#	C = fOne$rsCov
#	#
#	i = i+1
#}
#
##
#aS = fOne$aS
#a0 = fOne$a0
#M  = fOne$M
#kappa = fOne$kappa
#WW = fOne$WW
#ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
##
#B0 = 10000
#
##
#outlFV = getlFV(fOne)
#outlV = getlV(fOne)
#
##
#Dall = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
#D = Dall[Dall$lFV>0 & Dall$lB0V>0,]
#D = D[complete.cases(D),]
##who = (1/(D$xiHat+2)+0.001)>D$zetaHat & D$zetaHat>(1/(D$xiHat+2)-0.001) 
##D = D[!who,]
##
#D$lFV = D$lFV*fv #[D$lF>-3] = D$lFV[D$lF>-3]*10
##D = D[c(rep(T, 1), F),]
##D = D[seq(1, nrow(D), 2),]
##D = Dall[Dall$lF<4,]
##plot(D[,1], D[,2], pch=20)
##points(D[,3], D[,4])
##plot(D[,3], D[,4])
##points(D[,5], D[,5], pch=20)
#
##
##GP INTERPOLATION
##
#
##First fit the lF Model and make lF related Predictions
#
##pick a polynomial mean function
#lFy = D$lF
#lFV = diag(D$lFV)
#lFX = cbind(1, D$xiSeed, D$zetaSeed)
#
##
#xAug = seq(0.5, 4, 0.1) #0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
#aug = cbind(rep(1, length(xAug)), xAug, getZetaBH(xAug, M, kappa, WW, aS, a0)) #1/(xAug+2))
##xAug = numeric(0)
##aug = numeric(0)
#lFX = rbind(lFX, aug)
#lFy = c(lFy, log(xAug*M))
#lFVm = mean(D$lFV, na.rm=T)
#lFV = diag(c(D$lFV, rep(lFVm, length(xAug))))
#lFV[is.na(lFV)] = lFVm
#
##
#lFaxes = lFX[,2:3]
#registerData(lFy, lFX, lFaxes, lFV)
#par = c(l1=0.5, l2=0.5, s2=0.5)
#lFFit = gpMAP(par, hessian=F, psiSample=F, lower=c(rep(eps(), 3))) #lower=c(0.3, rep(eps(), 2)))
#writeLines("lF Fit:")
#print(lFFit)
#
##lF prediction
#zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.005) #length.out=) 	#seq(0.15, 0.35, 0.001)  #rev(seq(0.1, 0.80, 0.01)) #
#xiStar   = seq(min(D$xiSeed), max(D$xiSeed), 0.01) #length.out=)  	#seq(1, 3.5, 0.005)
#lFXStar = cbind(1, expand.grid(xiStar, zetaStar))
#mask = sapply(1:nrow(lFXStar), function(i){
#		#
#		win = 0.2#0.3
#		#
#		xi = lFXStar[i,2]
#		zeta = lFXStar[i,3]
#		#step = win*10/0.1
#		bot = c()
#		top = c()
#		for(x in seq(max(xi-win, xiBot), min(xi+win, xiTop), length.out=30)){
#			#
#			isWin = x<=(D$xiBin+win) & x>=(D$xiBin-win)
#			bot = c(bot, min(D$zetaBin[isWin]) )
#			top = c(top, max(D$zetaBin[isWin]) )
#		}
#		bot = mean(bot)
#		top = mean(top)
#		##
#		#bot = min(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#		#top = max(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#		#
#		return( zeta>bot & zeta<top & zeta>zetaBot & zeta<zetaTop )
#		#return(T)
#	}
#)
#lFPred = gpPredict(lFXStar, lFXStar[,2:3], lFFit)
#lFPred[!mask] = NA
#
##bias
#xiHat = exp(lFPred)/M
#xiBias = sweep(xiHat, 1, xiStar)
#xiBias[is.infinite(xiBias)] = NA
##zbh = BBar(exp(lFPred), M, kappa, ww, WW, fit$alpha, fit$beta, fit$gamma)/BBar(0, M, kappa, ww, WW, fit$alpha, fit$beta, fit$gamma) 
#zbh = matrix(getZetaBH(exp(lFPred)/M, M, kappa, WW, aS, a0), nrow(xiBias), ncol(xiBias))
#zetaBias = sweep(zbh, 2, zetaStar)
##zetaBias = sweep(matrix(0.5, nrow(xiHat), ncol(xiHat)), 2, zetaStar)
#zetaBias[!mask] = NA
#
##
#eucBias = mcmapply(function(xiHat, xi, zeta){
#                myDist(xiHat, xi, zeta)
#        }, xiHat, lFXStar[,2], lFXStar[,3], mc.cores=6 #detectCores()
#)
#eucBias = matrix(eucBias, nrow=length(xiStar), ncol=length(zetaStar))
#
##Now fit the lK Model and make lK related Predictions
#
##
#lKy = D$lB0
#lKV = diag(c(D$lB0V))
#lKX = cbind(1, D$xiSeed, D$zetaSeed)
#
###
##xAug = seq(0.5, 4, 0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
##aug = cbind(rep(1, length(xAug)), xAug, getZetaBH(xAug, M, kappa, WW, aS, a0)) #1/(xAug+2))
###xAug = seq(0.75, 4, 0.5) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
###aug = cbind(rep(1, length(xAug)), xAug, 1/(xAug+2))
##lKX = rbind(lKX, aug)
##lKy = c(lKy, log(B0))
##lKVm = mean(D$lB0V, na.rm=T)
##lKV = diag(c(D$lB0V, rep(lKVm, length(xAug))))
#
##
#lKaxes = lKX[,2:3]
##
#registerData(lKy, lKX, lKaxes, lKV)
#par = c(l1=0.5, l2=0.5, s2=0.5)
#lKFit = gpMAP(par, hessian=F, psiSample=F, lower=c(eps(), rep(eps(), 2)))
#writeLines("lK Fit:")
#print(lKFit)
#
##lK prediction
#lKXStar = cbind(1, expand.grid(xiStar, zetaStar))
#lKPred = gpPredict(lKXStar, lKXStar[,2:3], lKFit)
#lKPred[!mask] = NA
#
##biomass bias
#kBias = exp(lKPred)-B0
##zetaHat*B0Hat
#bMSYHat = zbh*exp(lKPred) #exp(lKPred)/(xiHat+2)
#bMSYBias = sweep(bMSYHat, 2, zetaStar*B0)
#
##MSY bias
#bMSYStar = zetaStar*B0
#fMSYStar = xiStar*M
###SRR(Bmsy, a, b, g)
##alpHat = getAlpha(-1, exp(lFPred), M)
##betHat = getBeta(alpHat, -1, M, exp(lKPred))
##msyHat = SRR(bMSYHat, alpHat, betHat, -1)
##msyTru = matrix(NA, nrow=length(xiStar), ncol=length(zetaStar))
###for(j in 1:length(zetaStar)){
###	#
###	for(i in 1:length(xiStar)){
###		#
###		par = getPar(xiStar[i], zetaStar[j], M)
###		msyTru[i, j] = SRR(bMSYStar[j], log(par[1]), log(P0), par[2])
###	}
###}
###msyBias = msyHat-msyTru
#
##
##PLOT
##
#
##F* bias
#
##
#freq = c(T,F,F,F,F,F)
#
##cluster
#
##
#bh = function(x){getZetaBH(x, M, kappa, WW, aS, a0)}
#
##
#eg = expand.grid(xiStar, zetaStar)
#ms = apply(eg, 1, function(r){ stats::optimize(function(x){ sqrt((r[1]-x)^2 + (r[2]-bh(x))^2) }, c(0,4))$objective })
#ms = matrix(ms, nrow=length(xiStar), ncol=length(zetaStar))
#
##
#cols = brewer.pal(5, 'Set1')
#
##
#library(GauPro)
#gp = GauPro(D$xiHat, D$xiSeed)
#gpThresh = ga( type='real-valued',
#        fitness=function(x){ -gp$predict(x) },
#        lower=0, upper=1, optim=T
#)@solution[1] 
##optim(0.3, gp$predict, method="L-BFGS-B")$par  #stats::optimize(gp$predict, c(0,3.5))$minimum  #0.4175116 
#gpse = function(x) gp$predict(x)+2*gp$predict(x, se=T)$se
##stats::optimize(gpse, c(0,3.5))$minimum
#metaMean = na.omit(rowMeans(xiHat, na.rm=T))[1]
#metaMedian = na.omit(rowMedians(xiHat, na.rm=T))[1]
##
#lFXSSmall = lFXStar[lFXStar[,2]==min(lFXStar[,2]),]
#seTry = sqrt(gpPredictVar(lFXSSmall, lFXSSmall[,2:3], tg=0, lFFit))
#xiHatLower = exp(lFPred[1,]-1.96*seTry)/M
#metaLower = mean(xiHatLower, na.rm=T)
#metaLowerMedian = median(xiHatLower, na.rm=T)
##metaLower = (exp(lFPred[1,]-1.96*seTry)/M)()
##var = SIGStar%*%strongInv(SIG)%*%
##maternCor(anisoNorm(axeStar, axeStar, gpFit$psi['l1'], gpFit$psi['l2'], gpFit$psi['th']), gpFit$psi['nu'])*gpFit$psi['s2'] -
##       ( maternCor(anisoNorm(axeStar, axes, gpFit$psi['l1'], gpFit$psi['l2'], gpFit$psi['th']), gpFit$psi['nu'])*gpFit$psi['s2'] ) %*%
##        strongInv( maternCor(anisoNorm(axes, axes, gpFit$psi['l1'], gpFit$psi['l2'], gpFit$psi['th']), gpFit$psi['nu'])*gpFit$psi['s2'] ) %*%
##       ( maternCor(anisoNorm(axes, axeStar, gpFit$psi['l1'], gpFit$psi['l2'], gpFit$psi['th']), gpFit$psi['nu'])*gpFit$psi['s2'] ) + 
##       Tg
#
##metaSE = sqrt(na.omit(rowVars(xiHat, na.rm=T))[1])
##metaThresh = metaMean-2*metaSE
##metaThresh = na.omit(rowMeans(xiHat, na.rm=T))[1]
##zbih
#
##
#zetaThresh0.5 = bh(0.5)
#diffs0.5 = abs(zetaStar-zetaThresh0.5)
#zetaThresh0.25 = bh(0.25)
#diffs = abs(zetaStar-zetaThresh0.25)
##metaLowerZetaThresh = xiHatLower[diffs==min(diffs)] #xiiHatLower[ zetaStar[diffs==min(diffs)] ]
##metaLowerZeta1SEThresh = (exp(lFPred[1,]-1*seTry)/M)[diffs==min(diffs)]
#metaLowerZetaThresh0.5 = (exp(lFPred[1,]-1.96*seTry)/M)[diffs0.5==min(diffs0.5)]
#metaLowerZeta2SEThresh0.5 = (exp(lFPred[1,]-2*seTry)/M)[diffs0.5==min(diffs0.5)]
#metaLowerZeta1SEThresh0.5 = (exp(lFPred[1,]-1*seTry)/M)[diffs0.5==min(diffs0.5)]
#metaLowerZeta0SEThresh0.5 = (exp(lFPred[1,]-0*seTry)/M)[diffs0.5==min(diffs0.5)]
#
##
#lFStar = log(xiStar*M)
#lFStar = apply(t(lFStar), 2, rep, ncol(lFPred))
##
#m = (lFPred-t(lFStar))/t(lFStar)
#var = gpPredictVar(lFXStar, lFXStar[,2:3], tg=0, lFFit)
##se = sqrt(var/t(lFStar^2))
###
##metaRel0SEThresh = exp(m)/M
##metaRel1SEThresh = exp(m-1*se)/M
##metaRel2SEThresh = exp(m-2*se)/M
#
#
##
###
##png(sprintf("metaLowerZetaLineDD%s.png", mod))
##image(xiStar, zetaStar, xiHat<metaLowerZeta0SEThresh0.5 & xiHat>metaLowerZeta2SEThresh0.5, #-0.09,
##        col = c("white","black"), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}/M$"),
##        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
##        main = TeX("BH Inference Failure"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##curve(bh(x), from=0, to=4, lwd=3, add=T)
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##dev.off()
##
####
###png(sprintf("metaLowerZeta1SEThreshDD%s.png", mod))
###image(xiStar, zetaStar, xiHat>metaLowerZeta1SEThresh, #zbh<=zetaThresh, #0.445,
###        col = cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
###        xlab = TeX("$F_{MSY}/M$"),
###        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
###        main = TeX("BH Inference Failure"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
###        ylim = c(zetaBot, zetaTop),
###        xlim = c(xiBot, xiTop),
###        cex.lab = 1.5,
###        cex.main= 1.5
###)
###curve(bh(x), from=0, to=4, lwd=3, add=T)
###points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
###dev.off()
##
###
##png(sprintf("metaLowerZeta1SEThresh0.5DD%s.png", mod))
##image(xiStar, zetaStar, xiHat>metaLowerZeta1SEThresh0.5, #zbh<=zetaThresh, #0.445,
##        col = cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}/M$"),
##        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
##        main = TeX("BH Inference Failure"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##curve(bh(x), from=0, to=4, lwd=3, add=T)
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##dev.off()
##
### 
##png(sprintf('pairs%s.png', mod))
##pairs(D[,c('xiSeed', 'zetaSeed', 'xiHat', 'zetaHat')])
##dev.off()
#
##
#save.image(sprintf("%s.RData", mod))

#
#COMBINE
#

##
#load('FlatT45N150A0-1AS2K0.1N84Edge.RData')
#
###c('red', 'purple', 'blue')
#cols = cols[c(1,4,2)]
##
##png(sprintf("metaLowerZetaLineDD%s.png", mod))
#chop = 1:25 #((ncol(xiHat)-30):ncol(xiHat))
#m = matrix(mask, nrow(xiHat), ncol(xiHat))
#m[,chop] = F
##
##(xiHat<metaLowerZeta0SEThresh0.5 & xiHat>metaLowerZeta2SEThresh0.5 & m)
##
#heavySide = (xiHat<metaLowerZeta0SEThresh0.5) 
#lineMat = edge(heavySide)
#bound = ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y
##
#heavySide = (xiHat<metaLowerZeta2SEThresh0.5) 
#lineMat = edge(heavySide)
#bound = c(bound, rev( ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y ))
##
#heavySide = (xiHat<metaLowerZeta1SEThresh0.5) # & xiHat>metaLowerZeta2SEThresh0.5 & m) #xiHat>metaLowerZeta2SEThresh0.5 & m
#lineMat = edge(heavySide)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##
#png(sprintf("metaLowerZetaLinesDD%s.png", mod))
#plot( ksmooth(xiStar, ksmooth(xiStar, line)$y), 
#	type = 'l', 
#	lwd  = 3, 
#	col  = cols[3], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon( c(xiStar, rev(xiStar)), bound, col=makeTransparent(cols[3]), border=NA)
#
##image(xiStar, zetaStar, line, #(xiHat<metaLowerZeta0SEThresh0.5 & xiHat>metaLowerZeta2SEThresh0.5 & m), #-0.09,
##        col = c(NA,cols[3]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}/M$"),
##        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
##        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
###curve(bh(x), from=0, to=4, lwd=3, add=T, col=cols[3])
###points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#
##
#load('FlatT45N150A0-1AS4K0.2N56.RData')
##
#heavySide = (xiHat<metaLowerZeta0SEThresh0.5)
#lineMat = edge(heavySide)
#bound = ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y
##
#heavySide = (xiHat<metaLowerZeta2SEThresh0.5)
#lineMat = edge(heavySide)
#bound = c(bound, rev( ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y ))
##
#heavySide = (xiHat<metaLowerZeta1SEThresh0.5) # & xiHat>metaLowerZeta2SEThresh0.5 & m) #xiHat>metaLowerZeta2SEThresh0.5 & m
#lineMat = edge(heavySide)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##png(sprintf("metaLowerZetaLineDD%s.png", mod))
##image(xiStar, zetaStar, line, #xiHat<metaLowerZeta0SEThresh0.5 & xiHat>metaLowerZeta2SEThresh0.5, #-0.09,
##        col = c(NA,cols[2]), #makeTransparent(c(NA,cols[2])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
#cols = cols[c(1,4,2)]
#lines( ksmooth(xiStar, ksmooth(xiStar, line)$y),
#        type = 'l',
#        lwd  = 3,
#        col  = cols[2],
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	#add=T
#)
#polygon( c(xiStar, rev(xiStar)), bound, col=makeTransparent(cols[2]), border=NA)
#
##curve(bh(x), from=0, to=4, lwd=3, add=T, col=cols[2])
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##dev.off()
#
##
#load('FlatT45N150A0-1AS0.1K10N56.RData')
##
#heavySide = (xiHat<metaLowerZeta0SEThresh0.5+0.04)
#lineMat = edge(heavySide)
#bound = ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y
##
#heavySide = (xiHat<metaLowerZeta2SEThresh0.5-0.04)
#lineMat = edge(heavySide)
#bound = c(bound, rev( ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y ))
##
#heavySide = (xiHat<metaLowerZeta1SEThresh0.5) # & xiHat>metaLowerZeta2SEThresh0.5 & m) #xiHat>metaLowerZeta2SEThresh0.5 & m
#lineMat = edge(heavySide)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##
##image(xiStar, zetaStar, line, #xiHat<metaLowerZeta0SEThresh0.5 & xiHat>metaLowerZeta2SEThresh0.5-0.09,
##        col = c(NA,cols[1]), #makeTransparent(c(NA,cols[1])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
#cols = cols[c(1,4,2)]
#lines( ksmooth(xiStar, ksmooth(xiStar, line)$y),
#        type = 'l',
#        lwd  = 3,
#        col  = cols[1],
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	#add=T
#)
#polygon( c(xiStar, rev(xiStar)), bound, col=makeTransparent(cols[1]), border=NA)
#
##curve(bh(x), from=0, to=4, lwd=3, add=T, col=cols[1])
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##dev.off()
#
##
#legend('topright', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, WW, a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, WW, a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, WW, a0),2)))), fill=cols[c(3,2,1)])
#
#dev.off()

##
#howBad = 0.1
#
##
#load('FlatT45N150A0-1AS2K0.1N84Edge.RData')
#
###
##heavySide = (xiHat<metaLowerZeta0SEThresh0.5) 
##lineMat = edge(heavySide)
##bound = ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y
###
##heavySide = (xiHat<metaLowerZeta2SEThresh0.5) 
##lineMat = edge(heavySide)
##bound = c(bound, rev( ksmooth(xiStar, ksmooth(xiStar,apply(lineMat, 1, function(x)zetaStar[x]))$y)$y ))
###
##heavySide = (xiHat<metaLowerZeta1SEThresh0.5) # & xiHat>metaLowerZeta2SEThresh0.5 & m) #xiHat>metaLowerZeta2SEThresh0.5 & m
##lineMat = edge(heavySide)
##line = apply(lineMat, 1, function(x)zetaStar[x])
###
##png(sprintf("metaLowerZetaLinesDD%s.png", mod))
##plot( ksmooth(xiStar, ksmooth(xiStar, line)$y), 
##	type = 'l', 
##	lwd  = 3, 
##	col  = cols[3], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}/M$"),
##        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
##        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##polygon( c(xiStar, rev(xiStar)), bound, col=makeTransparent(cols[3]), border=NA)
#
##
#image(xiStar, zetaStar, m-2*se>howBad,
#        col = c(NA,cols[3]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#
##
#load('FlatT45N150A0-1AS4K0.2N56.RData')
#image(xiStar, zetaStar, (m-2*se)>howBad,
#        col = c(NA,cols[2]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	add = T
#)
#
##
#load('FlatT45N150A0-1AS0.1K10N56.RData')
##
#var[var<0]=eps()
#var = fv*var
#se = sqrt(var/t(lFStar^2))
##
#image(xiStar, zetaStar, (m-2*se)>howBad,
#        col = c(NA,cols[1]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	add = T
#)
#
###
##load('FlatT45N150A0-1AS1K0.5N56.RData')
###
##image(xiStar, zetaStar, (m-2*se)>howBad,
##        col = c(NA,cols[1]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
##        xlab = TeX("$F_{MSY}/M$"),
##        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
##        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5,
##	add = T
##)
#
#
##
#howBad = 0.5
#
#dev.new()
##
#load('FlatT45N150A0-1AS2K0.1N84Edge.RData')
#image(xiStar, zetaStar, m-2*se>howBad,
#        col = c(NA,cols[3]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eucCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#
##
#load('FlatT45N150A0-1AS4K0.2N56.RData')
#image(xiStar, zetaStar, (m-2*se)>howBad,
#        col = c(NA,cols[2]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	add = T
#)
#
##
#load('FlatT45N150A0-1AS0.1K10N56.RData')
#var[var<0]=eps()
#var = fv*var
#se = sqrt(var/t(lFStar^2))
#image(xiStar, zetaStar, (m-2*se)>howBad,
#        col = c(NA,cols[1]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	add = T
#)
#
###
##load('FlatT45N150A0-1AS1K0.5N56.RData')
###
##image(xiStar, zetaStar, (m-2*se)>howBad,
##        col = c(NA,cols[1]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
##        xlab = TeX("$F_{MSY}/M$"),
##        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
##        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5,
##	add = T
##)

##
#lFStar = log(xiStar*M)
#lFStar = apply(t(lFStar), 2, rep, ncol(lFPred))
##
#m = (lFPred-t(lFStar))/t(lFStar)
#var = gpPredictVar(lFXStar, lFXStar[,2:3], tg=0, lFFit)
#se = sqrt(var/t(lFStar^2))
##
#metaRel0SEThresh = exp(m)/M
#metaRel1SEThresh = exp(m-1*se)/M
#metaRel2SEThresh = exp(m-2*se)/M

#
howBad = 0.6

##
#png(sprintf("relErrorLinesDD%s.png", howBad))
##
#load('FlatT45N150A0-1AS2K0.1N84Edge.RData')
#cols = cols[c(1,4,2)]
##
#var[var<0]=eps()
#var = var/fv#/fv
#se = sqrt(var)
##
#heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
#lineMat = edge(heavySide) #hsToLine(HS)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##
#plot( ksmooth(xiStar, ksmooth(xiStar, line)$y), 
#	type = 'l', 
#       	lwd  = 3, 
#       	col  = cols[3], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#
##
#load('FlatT45N150A0-1AS1K0.5N56.RData')
#cols = cols[c(1,4,2)]
#
##
#var[var<0]=eps()
#var = var/fv#/fv
#se = sqrt(var)
##
#heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
#lineMat = edge(heavySide) #hsToLine(HS)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##
#lines( ksmooth(xiStar, ksmooth(xiStar, line)$y), 
#	type = 'l', 
#       	lwd  = 3, 
#       	col  = cols[2], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#
##
#load('FlatT45N150A0-1AS0.1K10N56.RData')
#cols = cols[c(1,4,2)]
#
##
#var[var<0]=eps()
#var = var/fv#/fv
#se = sqrt(var)
##
#heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
#lineMat = edge(heavySide) #hsToLine(HS)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##
#lines( ksmooth(xiStar, ksmooth(xiStar, line)$y), 
#	type = 'l', 
#       	lwd  = 3, 
#       	col  = cols[1], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#
###
##legend('topright', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, WW, a0), 2))), TeX(sprintf("$a_s=1$,    $\\kappa=0.5$, $w(a_s)\\approx %s$", round(vbGrow(1, 0.5, WW, a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, WW, a0),2)))), fill=cols[c(3,2,1)])
#
#dev.off()

##
#load('FlatT45N150A0-1AS4K0.2N56.RData')
#cols = cols[c(1,4,2)]
##
#var[var<0]=eps()
#var = var/fv/fv
#se = sqrt(var)
##
#heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
#lineMat = edge(heavySide) #hsToLine(HS)
#line = apply(lineMat, 1, function(x)zetaStar[x])
##
#lines( ksmooth(xiStar, ksmooth(xiStar, line)$y), 
#	type = 'l', 
#       	lwd  = 3, 
#       	col  = cols[1], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)

#
png(sprintf("relErrorImagesBHDD%s.png", howBad))

#
load('FlatT45N150A0-1AS2K0.1N84Edge.RData')
cols = cols[c(1,4,2)]

#
var[var<0]=eps()
var = var/fv #var/fv/fv
se = sqrt(var)
#
heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
lineMat = edgeBH(heavySide) #hsToLine(HS)
line = apply(lineMat, 1, function(x)zetaStar[x])
#
image( xiStar, zetaStar, lineMat, 
	#type = 'l', 
       	lwd  = 3, 
       	col  = c(NA,cols[3]), #cols[3], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        cex.lab = 1.5,
        cex.main= 1.5
)

#
load('FlatT45N150A0-1AS1K0.5N56.RData')
cols = cols[c(1,4,2)]

#
var[var<0]=eps()
var = var/fv #var/fv/fv
se = sqrt(var)
#
heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
lineMat = edgeBH(heavySide) #hsToLine(HS)
line = apply(lineMat, 1, function(x)zetaStar[x])
#
image( xiStar, zetaStar, lineMat,
	#type = 'l', 
       	lwd  = 3, 
       	col  = c(NA,cols[2]), #cols[2], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        cex.lab = 1.5,
        cex.main= 1.5,
	add = T
)

#
load('FlatT45N150A0-1AS0.1K10N56.RData')
cols = cols[c(1,4,2)]

#
var[var<0]=eps()
var = var/fv#var/fv/fv
se = sqrt(var)
#
heavySide = qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar))
lineMat = edgeBH(heavySide) #hsToLine(HS)
line = apply(lineMat, 1, function(x)zetaStar[x])
#
image( xiStar, zetaStar, lineMat,
	#type = 'l', 
       	lwd  = 3, 
       	col  = c(NA,cols[1]), #cols[1], #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), #col  = adjustcolor(eu
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointl
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        cex.lab = 1.5,
        cex.main= 1.5,
	add = T
)

##
#legend('topright', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, WW, a0), 2))), TeX(sprintf("$a_s=1$,    $\\kappa=0.5$, $w(a_s)\\approx %s$", round(vbGrow(1, 0.5, WW, a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, WW, a0),2)))), fill=cols[c(3,2,1)])

dev.off()










#
#JUNK
#

##
#howBad = 0.5 #0.5 #0.1
#
#dev.new()
#load('FlatT45N150A0-1AS4K0.2N56.RData')
#load('FlatT45N150A0-1AS2K0.1N84Edge.RData')
#load('FlatT45N150A0-1AS0.1K10N56.RData')
#load('FlatT45N150A0-1AS1K0.5N56.RData')

#
##
#FPred = exp(lFPred+var/2)
#FVar = (exp(var)-1)*FPred^2
##
#xiPred = FPred/M
#xiVar = FVar/(M^2)
###
##m = (xiPred-t(exp(lFStar)/M))/t(exp(lFStar)/M)
##se = sqrt(xiVar/t(exp(lFStar)/M)^2)
#se = sqrt(var)
#
#var = fv*var
#se = sqrt(var/t(lFStar^2))



##
#image(xiStar, zetaStar, qlnorm(0.05, lFPred, se)<(1-howBad)*t(exp(lFStar)), #(m-2*se)>log(-howBad*t(exp(lFStar))+t(exp(lFStar))),
#        col = c(NA,cols[3]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)

##
#image(xiStar, zetaStar, (m-2*se)>log(-howBad*t(exp(lFStar))+t(exp(lFStar))),
#        col = c(NA,cols[3]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)

##
#howBad = 0.2
##
#image(xiStar, zetaStar, qlnorm(0.025, lFPred, sqrt(var))>(-howBad*t(exp(lFStar))+t(exp(lFStar))), #qlnorm(0.025, lFPred, sqrt(var))>howBad, #(m+2*se)<howBad,
#        col = c(NA,cols[2]), #makeTransparent(c(NA,cols[3])), #cols[1:2], #('red', 'green'), 
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("BH Inference Failure Thresholds"), #"Bias Direction for ($F_{MSY}/M$, B_{
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5,
#	add=T
#)






