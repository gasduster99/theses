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
#DATA
#

#
#mod = #"ExpT45N150K1" #"FlatT30N150A15K0.1" #"ExpT45N150A15" #"ExpT45N150A7.5" # #"ExpT45N150Wide"
#fail
#mod = "FlatT30N150A-1AS15K0.1"
#fail
#mod = "FlatT30N150A-0.5AS15K0.1" 
##worked w/ F help w/o K help
#mod = "ExpT45N150A-1AS2" 
##worked w/ F help w/o K help
#mod = "ExpT45N150A-1AS15K0.1" 
#fail
#mod = "ExpT45N150A-0.5AS15K0.1" 

## worked but no xi<1
#mod = "ExpT45N300AS1K1" 
#kinda worked
#mod = "ExpT45N300AS10K10"

##failed
#mod = "ExpT45N300AS10K0.1"
#WORKED
#mod = "ExpT45N300A0-1AS10K0.1"
#GOOD START: add a bit of refinement
#mod = "ExpT45N150A0-1AS2K0.1"
#failed
#mod = "ExpT45N300A0-1AS10K0.1N28"
##WORKED 
#mod = "ExpT45N300AS0.1K10"

##failed Fmsy overestimated
#mod = "FlatT45N300A0-1AS10K0.1"
#bumpy 
#mod = "FlatT45N300A0-1AS10K0.1N56"
#mod = "FlatT45N300A0-1AS10K0.1N84"
#failed		
#mod = "FlatT45N150A0-1AS3K0.1"
#GOOD START: add refinement
#mod = "FlatT45N150A0-1AS0.1K10"
#mod = "FlatT45N150A0-1AS0.1K10N28"
#mod = "FlatT45N150A0-1AS0.1K10N56"
#GOOD START: add refinment
#mod = "FlatT45N150A0-1AS2K0.1"
#mod = "FlatT45N150A0-1AS2K0.1N28"
#mod = "FlatT45N150A0-1AS2K0.1N56"
#mod = "FlatT45N150A0-1AS2K0.1N84"; fv=100
mod = "FlatT45N150A0-1AS2K0.1N84Edge"; fv=100

#
place = sprintf("./modsDD%s/", mod)

#
xiRes = 0.5
zetaTop = 0.6 #0.7
zetaBot = 0.2 #0.1
xiBot = 0.5
xiTop = 3.5

#
i = 1
C = NULL
while( is.null(C) ){
	#
	f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit_*.rda"))[i] )
	fOne = readRDS(f)
	C = fOne$rsCov
	#
	i = i+1
}

#
aS = fOne$aS
a0 = fOne$a0
M  = fOne$M
kappa = fOne$kappa
WW = fOne$WW
ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
#
B0 = 10000

#
outlFV = getlFV(fOne)
outlV = getlV(fOne)

#
Dall = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
D = Dall[Dall$lFV>0 & Dall$lB0V>0,]
D = D[complete.cases(D),]
#
D$lFV = D$lFV*fv #[D$lF>-3] = D$lFV[D$lF>-3]*10
#D = D[c(rep(T, 1), F),]
#D = D[seq(1, nrow(D), 2),]
#D = Dall[Dall$lF<4,]
#plot(D[,1], D[,2], pch=20)
#points(D[,3], D[,4])
#plot(D[,3], D[,4])
#points(D[,5], D[,5], pch=20)

#
#GP INTERPOLATION
#

#First fit the lF Model and make lF related Predictions

#pick a polynomial mean function
lFy = D$lF
lFV = diag(D$lFV)
lFX = cbind(1, D$xiSeed, D$zetaSeed)

#
xAug = seq(0.5, 4, 0.1) #0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
aug = cbind(rep(1, length(xAug)), xAug, getZetaBH(xAug, M, kappa, WW, aS, a0)) #1/(xAug+2))
#xAug = numeric(0)
#aug = numeric(0)
lFX = rbind(lFX, aug)
lFy = c(lFy, log(xAug*M))
lFVm = mean(D$lFV, na.rm=T)
lFV = diag(c(D$lFV, rep(lFVm, length(xAug))))
lFV[is.na(lFV)] = lFVm

#
lFaxes = lFX[,2:3]
registerData(lFy, lFX, lFaxes, lFV)
par = c(l1=0.5, l2=0.5, s2=0.5)
lFFit = gpMAP(par, hessian=F, psiSample=F, lower=c(rep(eps(), 3))) #lower=c(0.3, rep(eps(), 2)))
writeLines("lF Fit:")
print(lFFit)

#lF prediction
zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.005) #length.out=) 	#seq(0.15, 0.35, 0.001)  #rev(seq(0.1, 0.80, 0.01)) #
xiStar   = seq(min(D$xiSeed), max(D$xiSeed), 0.01) #length.out=)  	#seq(1, 3.5, 0.005)
lFXStar = cbind(1, expand.grid(xiStar, zetaStar))
mask = sapply(1:nrow(lFXStar), function(i){
		#
		win = 0.2#0.3
		#
		xi = lFXStar[i,2]
		zeta = lFXStar[i,3]
		#step = win*10/0.1
		bot = c()
		top = c()
		for(x in seq(max(xi-win, xiBot), min(xi+win, xiTop), length.out=30)){
			#
			isWin = x<=(D$xiBin+win) & x>=(D$xiBin-win)
			bot = c(bot, min(D$zetaBin[isWin]) )
			top = c(top, max(D$zetaBin[isWin]) )
		}
		bot = mean(bot)
		top = mean(top)
		##
		#bot = min(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
		#top = max(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
		#
		return( zeta>bot & zeta<top & zeta>zetaBot & zeta<zetaTop )
		#return(T)
	}
)
lFPred = gpPredict(lFXStar, lFXStar[,2:3], lFFit)
lFPred[!mask] = NA

#bias
xiHat = exp(lFPred)/M
xiBias = sweep(xiHat, 1, xiStar)
xiBias[is.infinite(xiBias)] = NA
#zbh = BBar(exp(lFPred), M, kappa, ww, WW, fit$alpha, fit$beta, fit$gamma)/BBar(0, M, kappa, ww, WW, fit$alpha, fit$beta, fit$gamma) 
zbh = matrix(getZetaBH(exp(lFPred)/M, M, kappa, WW, aS, a0), nrow(xiBias), ncol(xiBias))
zetaBias = sweep(zbh, 2, zetaStar)
#zetaBias = sweep(matrix(0.5, nrow(xiHat), ncol(xiHat)), 2, zetaStar)
zetaBias[!mask] = NA

#
eucBias = mcmapply(function(xiHat, xi, zeta){
                myDist(xiHat, xi, zeta)
        }, xiHat, lFXStar[,2], lFXStar[,3], mc.cores=6 #detectCores()
)
eucBias = matrix(eucBias, nrow=length(xiStar), ncol=length(zetaStar))

#Now fit the lK Model and make lK related Predictions

#
lKy = D$lB0
lKV = diag(c(D$lB0V))
lKX = cbind(1, D$xiSeed, D$zetaSeed)

##
#xAug = seq(0.5, 4, 0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
#aug = cbind(rep(1, length(xAug)), xAug, getZetaBH(xAug, M, kappa, WW, aS, a0)) #1/(xAug+2))
##xAug = seq(0.75, 4, 0.5) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
##aug = cbind(rep(1, length(xAug)), xAug, 1/(xAug+2))
#lKX = rbind(lKX, aug)
#lKy = c(lKy, log(B0))
#lKVm = mean(D$lB0V, na.rm=T)
#lKV = diag(c(D$lB0V, rep(lKVm, length(xAug))))

#
lKaxes = lKX[,2:3]
#
registerData(lKy, lKX, lKaxes, lKV)
par = c(l1=0.5, l2=0.5, s2=0.5)
lKFit = gpMAP(par, hessian=F, psiSample=F, lower=c(eps(), rep(eps(), 2)))
writeLines("lK Fit:")
print(lKFit)

#lK prediction
lKXStar = cbind(1, expand.grid(xiStar, zetaStar))
lKPred = gpPredict(lKXStar, lKXStar[,2:3], lKFit)
lKPred[!mask] = NA

#biomass bias
kBias = exp(lKPred)-B0
#zetaHat*B0Hat
bMSYHat = zbh*exp(lKPred) #exp(lKPred)/(xiHat+2)
bMSYBias = sweep(bMSYHat, 2, zetaStar*B0)

#MSY bias
bMSYStar = zetaStar*B0
fMSYStar = xiStar*M
##SRR(Bmsy, a, b, g)
#alpHat = getAlpha(-1, exp(lFPred), M)
#betHat = getBeta(alpHat, -1, M, exp(lKPred))
#msyHat = SRR(bMSYHat, alpHat, betHat, -1)
#msyTru = matrix(NA, nrow=length(xiStar), ncol=length(zetaStar))
##for(j in 1:length(zetaStar)){
##	#
##	for(i in 1:length(xiStar)){
##		#
##		par = getPar(xiStar[i], zetaStar[j], M)
##		msyTru[i, j] = SRR(bMSYStar[j], log(par[1]), log(P0), par[2])
##	}
##}
##msyBias = msyHat-msyTru

#
#PLOT
#

#F* bias

#
freq = c(T,F,F,F,F,F)

#
png(sprintf("fMSYBiasDD%s.png", mod))
nCols = 100
maxBias = abs(max(xiBias, na.rm=T))
minBias = abs(min(xiBias, na.rm=T))
posCols = hcl.colors(round(nCols*(maxBias/(maxBias+minBias))), "Reds 2", rev=T)
negCols = hcl.colors(round(nCols*(minBias/(maxBias+minBias))), "Blues 2", rev=F)
xCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, xiBias, 
	col  = adjustcolor(xCols, alpha.f=0.6),  
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'),
	main = TeX("Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop),
	xlim = c(xiBot, xiTop),
	cex.lab = 1.5,
        cex.main= 1.5
)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T) 
show = seq(1, length(xCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.1f", rev(seq(min(xiBias[xiMask, zetaMask], na.rm=T)*M, max(xiBias[xiMask, zetaMask], na.rm=T)*M, length.out=length(show)))),
#        fill = rev(xCols[show]), #colMap[c(1, 10, 20)], 
#        xpd = NA
#)
dev.off()

#F* relative bias

#
chaosThresh = 3.56995
rSurf = exp(lFPred)*2

#
png(sprintf("fMSYRelBiasDD%s.png", mod))
nCols = 21 #50*2
maxAbsXBias = 1 #max(abs(xiBias/xiStar), na.rm=T)
posCols = hcl.colors(nCols/2, "Reds 2", rev=T)
negCols = hcl.colors(nCols/2, "Blues 2", rev=F)
xCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, xiBias/xiStar, 
	col  = adjustcolor(xCols, alpha.f=0.6), 
        xlab = TeX("$F_{MSY}/M$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop),
	xlim = c(xiBot, xiTop),
	zlim = c(-maxAbsXBias, maxAbsXBias),
	cex.lab = 1.5,
        cex.main= 1.5
)
contour(xiStar[15:16], zetaStar, (xiBias/xiStar)[15:16,], 
	method = "simple", 
        levels = seq(-1, 1, 0.1),
	labcex = 1.025,
	lwd  = eps(),
	xlab = TeX("$F_{MSY}/M$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop), 
	xlim = c(xiBot, xiBot), 
	zlim = c(-maxAbsXBias, maxAbsXBias), 
	add  = T 
)
end = rev(length(xiStar)-c(15:16))
contour(xiStar[end], zetaStar, (xiBias/xiStar)[end,], 
	method = "simple", 
        levels = seq(-1, 1, 0.1),
	labcex = 1.025,
	lwd  = eps(),
	xlab = TeX("$F_{MSY}/M$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop), 
	xlim = c(xiBot, xiBot), 
	zlim = c(-maxAbsXBias, maxAbsXBias), 
	add  = T 
)
image(xiStar, zetaStar, xiBias/xiStar,
        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}/M$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(1, max(xiBias/xiStar, 1, na.rm=T)),
	add  = T
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
show = seq(1, length(xCols), length.out=nCols) #20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))), 
#        fill = rev(xCols[show]), 
#        xpd = NA
#)
points(D$xiSeed, D$zetaSeed)
dev.off()

#zeta bias

#
png(sprintf("zetaBiasDD%s.png", mod))
nCols = 100
maxBias = abs(max(zetaBias, na.rm=T))
minBias = abs(min(zetaBias, na.rm=T))
posCols = hcl.colors(round(nCols*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
negCols = hcl.colors(round(nCols*minBias/(maxBias+minBias)), "Blues 2", rev=F)
yCols = c(negCols, "#FFFFFF", posCols)
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, zetaBias,
	col  = adjustcolor(yCols, alpha.f=0.6), 
        xlab = TeX("$F_{MSY}/M$"), 
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Bias in Estimated $B_{MSY}/B_0$"), 
	ylim = c(zetaBot, zetaTop),
	xlim = c(xiBot, xiTop),
	cex.lab = 1.5,
        cex.main= 1.5
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
show = seq(1, length(yCols), length.out=20) 
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(zetaBias[xiMask, zetaMask], na.rm=T), max(zetaBias[xiMask, zetaMask], na.rm=T), length.out=length(show)))),
#        fill = rev(yCols[show]), 
#        xpd = NA
#)
points(D$xiSeed, D$zetaSeed)
dev.off()

#euc bias

#
png(sprintf("directionalBiasDD%s.png", mod))
eucCols = hcl.colors(41, "Reds 2", rev=T)
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, eucBias,
        col  = adjustcolor(eucCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}/M$"), 
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
	main = TeX("Bias Direction for ($F_{MSY}/M$, B_{MSY}/B_0) Jointly"),
	ylim = c(zetaBot, zetaTop),
	xlim = c(xiBot, xiTop),
	cex.lab = 1.5,
        cex.main= 1.5
)
#curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
thin = c(T,rep(F,125))#135))
quiver(
        lFXStar[w,2][thin], lFXStar[w,3][thin],
        xiBias[w][thin], zetaBias[w][thin],
        scale=0.025
)
show = seq(1, length(eucCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=length(show)))),
#        fill = rev(eucCols[show]), 
#        xpd = NA
#)
#points(D$xiSeed, D$zetaSeed)
dev.off()

#K bias

#
png(sprintf("kBiasDD%s.png", mod))
#
maxBias = abs(max(kBias, na.rm=T))
minBias = abs(min(kBias, na.rm=T))
nCol = 25
posCols = hcl.colors(round(nCol*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
negCols = hcl.colors(round(nCol*minBias/(maxBias+minBias)), "Blues 2", rev=F)
kCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, kBias, #xiStar[xiMask], zetaStar[zetaMask], xBias[xiMask, zetaMask],i
        col  = adjustcolor(kCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Bias in Estimated $B_0$"), 
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
	cex.lab = 1.5,
        cex.main= 1.5
)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
show = seq(1, length(kCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"),
#        sprintf("%1.0f", round(rev(seq(min(kBias[xiMask, zetaMask], na.rm=T)*M, max(kBias[xiMask, zetaMask], na.rm=T)*M, length.out=length(show)))/5)*5),
#        fill = rev(kCols[show]),
#        xpd = NA
#)
dev.off()

#K relative bias

#
png(sprintf("kRelBiasDD%s.png", mod))
#
nCol = 21 #15*2
maxAbsXBias = 1 #max(abs(kBias/P0), na.rm=T)
posCols = hcl.colors(nCol/2, "Reds 2", rev=T)
negCols = hcl.colors(nCol/2, "Blues 2", rev=F)
xCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, kBias/B0, 
        col  = adjustcolor(xCols, alpha.f=0.6), 
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $B_0$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(-maxAbsXBias, maxAbsXBias),
	cex.lab = 1.5,
        cex.main= 1.5
)
contour(xiStar[15:16], zetaStar, (kBias/B0)[15:16,], 
	method = "simple", 
        levels = seq(-1, 1, 0.1),
	labcex = 1.025,
	lwd  = eps(),
	xlab = TeX("$F_{MSY}/M$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop), 
	xlim = c(xiBot, xiBot), 
	zlim = c(-maxAbsXBias, maxAbsXBias), 
	add  = T 
)
end = rev(length(xiStar)-c(15:16))
contour(xiStar[end], zetaStar, (kBias/B0)[end,], 
	method = "simple", 
        levels = seq(-1, 1, 0.1),
	labcex = 1.025,
	lwd  = eps(),
	xlab = TeX("$F_{MSY}/M$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop), 
	xlim = c(xiBot, xiBot), 
	zlim = c(-maxAbsXBias, maxAbsXBias), 
	add  = T 
)
image(xiStar, zetaStar, kBias/B0,
        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}/M$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(1, max(kBias/B0, 1, na.rm=T)),
        add  = T
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
show = seq(1, length(xCols), length.out=nCol) #20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
#        fill = rev(xCols[show]), 
#        xpd = NA
#)
dev.off()

#B* bias

#
png(sprintf("bMSYBiasDD%s.png", mod))
#
maxBias = abs(max(bMSYBias, na.rm=T))
minBias = abs(min(bMSYBias, na.rm=T))
nCol = 25
posCols = hcl.colors(round(nCol*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
negCols = hcl.colors(round(nCol*minBias/(maxBias+minBias)), "Blues 2", rev=F)
bMSYCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, bMSYBias, 
        col  = adjustcolor(bMSYCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Bias in Estimated $B_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
	cex.lab = 1.5,
        cex.main= 1.5
)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
show = seq(1, length(bMSYCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"), 
#        sprintf("%1.0f", round(rev(seq(min(bMSYBias[xiMask, zetaMask], na.rm=T)*M, max(bMSYBias[xiMask, zetaMask], na.rm=T)*M, length.out=length(show)))/5)*5),
#        fill = rev(bMSYCols[show]),
#        xpd = NA
#)
dev.off()

#B* relative bias

#
png(sprintf("bMSYRelBiasDD%s.png", mod))
#
nCols = 21 #2*15
maxAbsXBias = 1 #max(abs(t(bMSYBias)/(P0*zetaStar)), na.rm=T)
posCols = hcl.colors(nCols/2, "Reds 2", rev=T)
negCols = hcl.colors(nCols/2, "Blues 2", rev=F)
xCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, t(t(bMSYBias)/(B0*zetaStar)),
        col  = adjustcolor(xCols, alpha.f=0.6), 
        xlab = TeX("$F_{MSY}/M$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $B_{MSY}$"), 
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop)
        ,zlim = c(-maxAbsXBias, maxAbsXBias),
	cex.lab = 1.5,
        cex.main= 1.5
)
contour(xiStar[15:16], zetaStar, (t(t(bMSYBias)/(B0*zetaStar)))[15:16,], 
	method = "simple", 
        levels = seq(-1, 1, 0.1),
	labcex = 1.025,
	lwd  = eps(),
	xlab = TeX("$F_{MSY}/M$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop), 
	xlim = c(xiBot, xiBot), 
	zlim = c(-maxAbsXBias, maxAbsXBias), 
	add  = T 
)
end = rev(length(xiStar)-c(15:16))
contour(xiStar[end], zetaStar, (t(t(bMSYBias)/(B0*zetaStar)))[end,], 
	method = "simple", 
        levels = seq(-1, 1, 0.1),
	labcex = 1.025,
	lwd  = eps(),
	xlab = TeX("$F_{MSY}/M$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'), 
	main = TeX("Relative Bias in Estimated $F_{MSY}$"), 
	ylim = c(zetaBot, zetaTop), 
	xlim = c(xiBot, xiBot), 
	zlim = c(-maxAbsXBias, maxAbsXBias), 
	add  = T 
)
image(xiStar, zetaStar, t(t(bMSYBias)/(B0*zetaStar)),
        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}/M$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(1, max(t(t(bMSYBias)/(B0*zetaStar)), 1, na.rm=T)),
        add  = T
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
show = seq(1, length(xCols), length.out=nCol) #20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
#        fill = rev(xCols[show]),
#        xpd = NA
#)
dev.off()

##MSY bias
#
##
#png(sprintf("msyBiasSchnute%s.png", mod))
##
#maxBias = abs(max(msyBias, na.rm=T))
#minBias = abs(min(msyBias, na.rm=T))
#nCol = 25
#posCols = hcl.colors(round(nCol*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
#negCols = hcl.colors(round(nCol*minBias/(maxBias+minBias)), "Blues 2", rev=F)
#msyCols = c(negCols, "#FFFFFF", posCols)
##
#xiMask = xiStar>xiBot & xiStar<xiTop
#zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
##
#par(mar=c(5, 4, 4, 5)+0.1)
#image(xiStar*M, zetaStar, msyBias,
#        col  = adjustcolor(msyCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Bias in Estimated MSY"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot*M, xiTop*M)
#)
##curve(0.5, from=0, to=12, lwd=3, add=T) #col=map2color(0, hcl.colors(41, "RdBu", rev=T)),
#abline(h=0.5, lwd=3)
#show = seq(1, length(msyCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"), 
#        sprintf("%1.0f", round(rev(seq(min(msyBias[xiMask, zetaMask], na.rm=T)*M, max(msyBias[xiMask, zetaMask], na.rm=T)*M, length.out=length(show)))/5)*5),
#        fill = rev(msyCols[show]),
#        xpd = NA
#)
#dev.off()
#
##MSY relative bias
#
##
#png(sprintf("msyRelBiasSchnute%s.png", mod))
##
#nCols = 2*15
#maxAbsXBias = max(abs(msyBias/msyTru), na.rm=T)
#posCols = hcl.colors(nCols/2, "Reds 2", rev=T)
#negCols = hcl.colors(nCols/2, "Blues 2", rev=F)
#xCols = c(negCols, "#FFFFFF", posCols)
##
#xiMask = xiStar>xiBot & xiStar<xiTop
#zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
##
#par(mar=c(5, 4, 4, 5)+0.1)
#image(xiStar*M, zetaStar, msyBias/msyTru,
#        col  = adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated MSY"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot*M, xiTop*M)
#        ,zlim = c(-maxAbsXBias, maxAbsXBias)
#)
##curve(0.5, from=0, to=12, lwd=3, add=T) #col=map2color(0, hcl.colors(41, "RdBu", rev=T)),
#abline(h=0.5, lwd=3)
#show = seq(1, length(xCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
#        fill = rev(xCols[show]),
#        xpd = NA
#)
#dev.off()

#catch

##
#png(sprintf("catch%s.png", mod))
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[5], zetaSims[5])
#dat = readRDS(fileDat)
#plot(dat$time, dat$catch, xlab="Time", ylab="F(t)/F*", ylim=c(0, 2), main="Catch=F* (F(t)/F*) B(t)", type="l", lwd=3)
#dev.off()

#design
png(sprintf("designLine%s.png", mod))
plot(D$xiSeed, D$zetaSeed, xlab=TeX("$F_{MSY}/M$"), ylab=TeX('$B_{MSY}/B_0$'), main="Design", cex.lab = 1.5, cex.main= 1.5, xlim=c(xiBot, xiTop), ylim=c(zetaBot, zetaTop))
curve(1/(x+2), from=0, to=4, lwd=3, add=T)
dotX = lFXStar[!mask,2][freq]
dotY = lFXStar[!mask,3][freq]
points(dotX[dotY>zetaBot & dotY<zetaTop], dotY[dotY>zetaBot & dotY<zetaTop], pch='.')
#plot(D$xiSeed, D$zetaSeed, xlab=TeX("$F_{MSY}/M$"), ylab=TeX('$B_{MSY}/B_0$'), main="Design", cex.lab = 1.5, cex.main= 1.5)
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
dev.off()

#design
png(sprintf("design%s.png", mod))
plot(D$xiSeed, D$zetaSeed, xlab=TeX("$F_{MSY}/M$"), ylab=TeX('$B_{MSY}/B_0$'), main="Design", cex.lab = 1.5, cex.main= 1.5, xlim=c(xiBot, xiTop), ylim=c(zetaBot, zetaTop))
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
dotX = lFXStar[!mask,2][freq]
dotY = lFXStar[!mask,3][freq]
points(dotX[dotY>zetaBot & dotY<zetaTop], dotY[dotY>zetaBot & dotY<zetaTop], pch='.')
dev.off()


#
#JUNK
#

##
#PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) } #((alpha/(M+
#
##
#getBeta = function(alpha, gamma, M, B0){
#        (1-(M/alpha)^gamma)/B0/gamma
#}
#b = Vectorize(getBeta, c("alpha", "gamma"))
#
##
#getAlpha = function(gamma, ff, M){
#        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
#}
#a = Vectorize(getAlpha, "gamma")
#
##
#getZeta = function(gamma, ff, M){
#        #fgfm = ff*gamma/(ff+M)
#        #fgfm/(1 + fgfm - (M/ff+M)^gamma)
#
#        #
#        (1-((M+ff)/getAlpha(gamma, ff, M))^gamma) / (1-(M/getAlpha(gamma, ff, M))^gamma)
#}
#z = Vectorize(getZeta, "gamma")
#
##
#FMsy = function(alpha, gamma, M){
#        #
#	if(gamma==-1){ return(sqrt(alpha*M)-M) }
#	#
#        FUpper = alpha-M
#        root = uniroot.all(function(ff){ a(gamma, ff, M)-alpha }, c(0, FUpper))
#
#        #root = uniroot.all(function(ff){ 1 - ff/gamma/(M+ff) - (alpha/(M+ff))^(-1/gamma) }, c(0, 
#        ##
#        #if(length(root)<1){ return(NA) }
#        #if(length(root)>1){ return(root[1]) }
#        ##
#        return(root)
#}












#
#OLD 
#

##residuals
#take = X[,3]!=0.5
#yRes = y[take]
#XRes = X[take,]
#gpPredRes = gpPredict(XRes, XRes[,2:3], gpFit, asMat=F)
#png(sprintf('pellaResiduals%s.png', mod))
#plot(XRes[,2], XRes[,3], col=map2color(yRes-gpPredRes, hcl.colors(10, "Blue-Red 3", rev=F)), pch=19)
#dev.off()

##
#png(sprintf("pellaRes%sHist2.png", mod))
#hist(yRes-gpPredRes)
#dev.off()
#
##
#sqLoss = function(y, X, yStar, yVar, k){
#        #y      : n-vector of observed values  
#	#X      : nx2 matrix matched with y
#	#yStar  : n-vector of predictive means at each observation
#        #xStar  : n-vector of predictive variance at each observation
#        #k      : G&G decision weighting factor
#        #
#        #value: return squared error posterior predictive loss  
#
#	#
#	out = sum(yVar + (k/(k+1))*(y-yStar)^2)
#	return(out)
#		
#        ##
#        #xiStar = X[,2]
#        #zetaStar = X[,3]
#        ##
#        #ppl = 0
#        #for(i in 1:length(y)){
#        #        ##
#        #        #where = which(xStar[,'x^1y^0']==lons[i] & xStar[,'x^0y^1']==lats[i])
#        #        #ys = yStar[where,]
#        #        ##
#        #        ppl = ppl + yVar[i] + (k/(k+1))*(y[i]-yStar[i])^2
#        #        #pp = rpois(MM, theta*exp(beta*x[i]))
#        #        #ppl = ppl + var(pp) + (k/(k+1))*(y[i]-mean(pp))^2
#        #}
#        ##
#        #return( ppl )
#}
#

###
##bub = 0.2
##D = Dall[Dall$xiInv<=max(xiSims)*(1+bub),]
##D = D[D$xiInv>=min(xiSims)*(1-bub),]
##D = D[D$zetaInv<=max(zetaSims)*(1+bub),]
##D = D[D$zetaInv>=min(zetaSims)*(1-bub),]
##D = D[D$lFV!=0 & D$xiBH<20,] #lalpha==0.04280697 is a numerical issue
#cut = 300#10
#D = D[D$xiInv<cut,]
#D = D[D$zetaInv>0,]
#D = D[!is.na(D$xiInv),]
###
#png(sprintf('pellaDat%s.png', mod))
#plot(D[,c("xiInv", "zetaInv")], 
#	ylim=c(min(D[,c("zetaInv", "zetaHat")]), max(zetaTop, max(D[,c("zetaInv", "zetaHat")]))), 
#	xlim=c(min(D[,c("xiInv", "xiHat")]), xiTop), #max(D[, c("xiInv")])), 
#	main='Schnute Tomlinson',
#	col=map2color(D$minDist, hcl.colors(60, "Zissou 1", rev=T)),
#	pch=20
#)
#points(D[,c("xiHat", "zetaHat")], col=map2color(D$minDist, hcl.colors(60, "Zissou 1", rev=T)))
#dev.off()

#
#OLD OLD
#

#X = cbind(D$xi, D$zeta)
#
##
#xiPred   = seq(0.7,3.5,length.out=50)	#seq(0.5, 3.5, 0.01) 		#	
#zetaPred = seq(0.2,0.75,length.out=50)	#seq(0.1, 0.75, 0.01)   	#
#predMesh = as.matrix(expand.grid(xiPred, zetaPred))
#colnames(predMesh) = c('xi', 'zeta')
##
#line = lm(cllhs~xi+zeta, data=D)
##gp = gpModel$new( S2=S2, X=X, Y=D$cllhs, obsV=D$cllhsV, #lm=line, 
##        B  = line$coeff,	#c(4.1894, -0.1452, -10.3066), 
##	v0 = 0.343979394716728,
##        v1 = 0.0270240223433852,
##        s2 = 2.93501281,
##        cr = -0.433574770896697
##)
###
##optOut0 = gp$fit(c('v0', 'v1', 's2', 'cr'),
##        lower   = c(eps(), eps(), eps(), -1),
##        upper   = c(Inf, Inf, Inf, 1),
##        cov     = T
##)
##gp$save('gpHotStart.rda')
#gp = readRDS('gpHotStart.rda')
#gp$printSelf()
##
#bnd = rep(Inf, length(line$coeff))
#optOut1 = gp$fit(c('B'),
#        lower   = c(-bnd),
#        upper   = c(bnd),
#        cov     = T
#)
#gp$save('gpHotterStart.rda')
#gp$printSelf()
##
#bnd = rep(Inf, length(line$coeff))
#optOut2 = gp$fit(c('v0', 'v1', 's2', 'cr', 'B'),
#        lower   = c(eps(), eps(), eps(), -1, -bnd),
#        upper   = c(Inf, Inf, Inf, 1, bnd),
#	#gaBoost = list(run=30, popSize=detectCores(), parallel=detectCores()),
#        cov     = T
#)
#gp$save('gpBigOpt.rda')
#gp$printSelf()
#
##
#yy = gp$predictMean(predMesh)
##SS = diag(gp$predictVar(predMesh))
###
##xS = exp(2*yy+SS)*(exp(SS)-1)/M^2
##K = 10^6
##yS = mcmapply(function(yy, SS){
##		logitM(2, log(2*M)-yy, SS, K=K) - logitM(1, log(2*M)-yy, SS, K=K)^2
##	}, yy, SS, mc.cores=8
##)
##xS = matrix(xS, nrow=length(xiPred))
##yS = matrix(yS, nrow=length(xiPred))
#
##
#xDist = matrix((exp(yy)/M)-xiPred, nrow=length(xiPred)) #, byrow=T)
#yDist = matrix((1/(exp(yy)/M+2))-zetaPred, nrow=length(xiPred))#, byrow=T)
##
#bigD = mcmapply(function(xiHat, xi, zeta){
#                dist(xiHat, xi, zeta)
#        }, exp(yy)/M, xiPred, zetaPred, mc.cores=detectCores()
#)
#
##
##OUTPUT
##
#
##
#m = 2.5
#
##
#pdf('gpXiBiasFineRB.pdf')
##dev.new()
#filled.contour(xiPred, zetaPred, xDist, 
#	zlim=c(-m, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#		#
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#	},
#	plot.title = title(
#		main = "Bias in Estimated Optimal Fishing Rate",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#pdf('gpZetaBiasFineRB.pdf')
##dev.new()
#filled.contour(xiPred, zetaPred, yDist, 
#	zlim=c(-0.55,0.55),
#        ylab = expression(frac(F[msy]/F[0])),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')	
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#	},
#	plot.title = title(
#		main = "Bias in Estimated Optimal Biomass",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3"),
#)
#dev.off()
#
##
#pdf('gpBiasArrowB.pdf')
##dev.new()
#datGrid = expand.grid(xiSims, zetaSims)
#colnames(datGrid) = c('xi', 'zeta')
#w = (predMesh[,'xi']>0.5 & predMesh[,'xi']<3.5 & predMesh[,'zeta']>0.2 & predMesh[,'zeta']<0.75) #paste0(predMesh[,'xi'], predMesh[,'zeta'])%in%paste0(datGrid[,'xi'], datGrid[,'zeta']) #
#thin = c(T,rep(F,3))
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred)),
#        zlim=c(0, m),
#        plot.axes = {
#                quiver(
#			predMesh[w,'xi'][thin], predMesh[w,'zeta'][thin],
#			xDist[w][thin], yDist[w][thin],
#			scale=0.05
#		) 
#                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#        },
#	plot.title = title(
#		main = "Directional Bias",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	#key.title = title(main = "Directional Bias"),
#        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
#)
#dev.off()


##
#source("gpClass0.0.1.r")



##X's are 2-vectors; L is 2x2; s2 is scalar
#S2 = function(X0, X1, s2, v0, v1, cr){
#	maxD = dbinorm(X0[,1], X0[,2], X0[,1], X0[,2], v0, v1, sqrt(v0*v1)*cr, log=T)
#	s2*mcmapply(function(x01, x02, m){
#		exp(dbinorm(X1[,1], X1[,2], x01, x02, v0, v1, sqrt(v0*v1)*cr, log=T)-m)
#	}, X0[,1], X0[,2], maxD, mc.cores=detectCores())	
#}
#
##
#logitM = function(n, mu, sig, K=10^5){
#	mean( inv.logit(qnorm((1:(K-1))/K, mu, sig))^n )
#}













#
#dev.new()
#filled.contour(xiPred, zetaPred, yDist/sqrt(yS), 
#	zlim=c(-0.6,0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)




##all parameters scalar
#Kr = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
#K = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
#Kmat = function(X, l){ 
#	#X	:a vector of predictors
#	#l	:a scalar length scale
#	
#	#
#	sapply(X, function(x){K(x,X,l)})
#}
##
#Kn = function(X, Y, lx, ly){ 
#	#X	:a vector of predictors
#	#Y	:a vector of predictors
#	#lx	:a scalar length scale for X
#	#ly	:a scalar length scale for Y
#	
#	#
#	sapply(Y, function(y){K(y,Y,ly)})*sapply(X, function(x){K(x,X,lx)})
#}
##
#Ku = function(X, Y, lx, ly){ 
#	#X	:a vector of predictors
#	#Y	:a vector of predictors
#	#lx	:a scalar length scale for X
#	#ly	:a scalar length scale for Y
#	
#	#
#	sapply(Y, function(y){K(y,Y,ly)})+sapply(X, function(x){K(x,X,lx)})
#}
##
#KnPred = function(Xs, Ys, X, Y, lx, ly){
#	sapply(Y, function(y){K(y,Ys,ly)})*sapply(X, function(x){K(x,Xs,lx)})
#}
#
##
#loglikeL = function(y, D, lx, ly, r, s2){	#function(y, D, L, sig2){ #
#	#
#	X = cbind(D$xi, D$zeta)
#	SIG2 = S2(X, X, lx, ly, r, s2) + diag(D$cllhsV)
#		
#	#
#	return( dmvnorm(y, sigma=SIG2, log=T) )
#}





#pp = predict( line )
#ppPred = predict(line, predMesh)


###
##optOut = c(mean(c(0.34395317, 0.02705277)), 2.93501281) 
##names(optOut) = c('lx', 'sig2')
#optOut = c(0.34395317, 0.02705277, -0.43195056, 2.93501281) 
#names(optOut) = c('lx', 'ly', 'cr', 'sig2')
##optOut = c(0.34395317, 0.02705277, 2.93501281)
##names(optOut) = c('lx', 'ly', 'sig2') 
##time: 340
#t = system.time({ 
#	optOut = optim(optOut, 
#		function(x){ 	
#			-loglikeL(D$cllhs-pp, D, x[1], x[2], x[3], x[4])
#		},
#		lower = c(10^-6, 10^-6, -1, 0),
#		upper = c(Inf, Inf, 1, Inf),
#		method = "L-BFGS-B",
#		hessian = T
#	)
#})
#writeLines(sprintf("Opt: %s",t[3]))
##IL = matrix(c(optOut$par[1], optOut$par[3], optOut$par[3], optOut$par[2]), nr=2)
#lx   = optOut$par['lx']
#ly   = optOut$par['ly']
#r    = optOut$par['cr']
#sig2 = optOut$par['sig2']
## 
##KInv = chol2inv(chol( sig2*Kn(D$xi,D$zeta,lx,ly) + diag(D$cllhsV) ))
#X = cbind(D$xi, D$zeta)
#KInv = chol2inv(chol( S2(X, X, lx, ly, r, sig2) + diag(D$cllhsV) ))
##KPX  = sig2*KnPred(predMesh$xi, predMesh$zeta, D$xi, D$zeta, lx, ly)
#KPX = S2(X, predMesh, lx, ly, r, sig2)
##KXP  = sig2*KnPred(D$xi, D$zeta, predMesh$xi, predMesh$zeta, lx, ly)
##KPP  = sig2*KnPred(predMesh$xi, predMesh$zeta, predMesh$xi, predMesh$zeta, lx, ly)
####time: 66
#ys = ppPred + KPX%*%KInv%*%(D$cllhs-pp) #cloglog(KPX%*%KInv%*%(D$cllhs-pp), inverse=T) #
###time: 470
##Ss = KPP-KPX%*%KInv%*%KXP
###writeLines(sprintf("Krig Var: %s",t[3]))

##
##CLLHS OUTPUT
##
#
##
#pdf('gpCllhsFineR.pdf')
#filled.contour(xiPred, zetaPred, matrix(ys, nrow=length(xiPred), byrow=T), 
#	zlim=c(-3.5, 1.5), #c(-3,0),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	}
#)
#dev.off()
#
##
#m = 3.5
#
##
#pdf('gpXiBiasFineR.pdf')
#xDist = matrix((exp(ys)/M)-xiPred, nrow=length(xiPred), byrow=T)
##xDist[xDist<(-1.2)] = -10
#filled.contour(xiPred, zetaPred, xDist, 
#	zlim=c(-m, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#pdf('gpZetaBiasFineR.pdf')
#yDist = matrix((1/(exp(ys)/M+2))-zetaPred, nrow=length(xiPred), byrow=T)
#filled.contour(xiPred, zetaPred, yDist, 
#	zlim=c(-0.6,0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#bigD = mapply(function(xiHat, xi, zeta){ 
#		dist(xiHat, xi, zeta) 
#	}, exp(ys)/M, xiPred, zetaPred
#) 
##
#pdf('gpBiasFineR.pdf')
##m = 3.5 #quantile(abs(bigD), 0.5)
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T), 
#	zlim=c(0, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blues", rev = TRUE)
#)
#dev.off()


###
##w = predMesh$X>0.5 & predMesh$X<3.5 & predMesh$Y>0.2 & predMesh$Y<0.75
###
##pdf('gpBiasArrow.pdf')
##filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T),
##        zlim=c(0, m),
##        plot.axes = {
##                quiver(predMesh$X[w[c(T,F)]], predMesh$Y[w[c(T,F)]], xDist[w[c(T,F)]], yDist[w[c(T,F)]])
##                #for(i in 1:length(predMesh$xi)){
##                #       arrows(predMesh$xi[i], predMesh$zeta[i], predMesh$xi[i] + scale * u, predMesh$zeta 
##                #}
##                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
##        },
##        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
##)
##dev.off()
#
#
#
#
#
##sapply(1:length(xiPred),
##	function(xi){sapply(1:length(zetaPred), 
##		function(zeta){ dist((exp(ys)/M)[], xiPred[], zetaPred[]) }
##	)
##}
#
###
##pdf('gpCllhsHessVar.pdf')
##filled.contour(xiPred, zetaPred, matrix(diag(Ss), nrow=length(xiPred), byrow=T), 
##	zlim=c(0, 0.6),
##	plot.axes = {
##		points(D$xi, D$zeta)
##		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
##	}
##)
##dev.off()
#
#
#
#
###
##p = scatterplot3d(D$xi, D$zeta, D$cllhs)
##plane = lm(cllhs~xi+zeta, data=D)
##p$plane3d(plane)
##pp = predict(plane)
##scatterplot3d(D$xi, D$zeta, D$cllhs-pp)
#
###
###DISTANCE OUTPUT
###
##
###
##plane = lm(log(distHat-distMin)~xi+zeta, data=D)
##pp = predict(plane)
##yy = log(D$distHat-D$distMin)-pp
###
##yys = KPX%*%KInv%*%(yy)
###
##pdf('gpLogBiasNoHess.pdf')
##filled.contour(xiPred, zetaPred, matrix(yys, nrow=length(xiPred), byrow=T),
##        #zlim=c(-50, 20),
##        plot.axes = {
##                points(D$xi, D$zeta)
##                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3) 
##                #points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
##                #points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
##        }
##)
##dev.off()

##
#SRR = function(P, lalpha, lbeta, gamma){
#        #
#        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
#}
#
##
#PBar = function(ff, alpha, beta, gamma, M){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }
#
##
#FMsy = function(alpha, gamma, M){
#        ##
#        #FUpper = P0 #alpha-M
#        ##root = uniroot.all(function(ff){  ((gamma-1)/alpha)^(gamma-1) - (M+ff)^(
#        #root = uniroot.all(function(ff){ 1 - ((M+ff)*((gamma-1)/alpha))^(1/(gamma
#        ##
#        #if(length(root)<1){ return(NA) }
#        #if(length(root)>1){ return(root[1]) }
#        #
#        root = alpha/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
#        return(root)
#}
#
##
#getPar = function(xi, zeta, M){
#        #
#        gamma = 1/zeta
#        odd = (1-zeta)/zeta
#        alpha = xi*M*odd*(1-zeta)^(-odd)
#        #
#        return(c(alpha, gamma))
#}
#
##
#getData = function(dir, xiSims, zetaSims){
#	#dir	: a directory containing data
#	#xiSims	: xis of simulated data
#	#zetaSims: zetas of simulated data
#	#
#	#Seed values are the initiating values to launch the inversion
#	#Inv values are the numerical inversions actually found for the 3 parameter curve at each seed value
#	#BH values are the MLE fits under the BH model.
#	
#	#
#	di = 1
#	D = data.frame(xiSeed=double(), zetaSeed=double(), xiInv=double(), zetaInv=double(), xiHat=double(), zetaHat=double(), minDist=double(), lF=double(), lFV=double(), lK=double(), lKV=double(), stringsAsFactors=F)
#	for(i in 1:length(zetaSims)){
#        	for(j in 1:length(xiSims)){
#                	#
#                	fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
#                	fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
#			
#			#
#			if( file.exists(fileDat) & file.exists(fileFit) ){
#				#
#				dat = readRDS(fileDat)
#				fit = readRDS(fileFit)	
#		
#				#
#				if(length(dat$time)!=length(dat$N) | length(fit$time)!=length(fit$N)){ next }
#						
#				#the inversion actually found
#				FInv = FMsy(dat$alpha, dat$gamma, M)
#				xiInv = FInv/M
#				zetaInv = PBar(FInv, dat$alpha, dat$beta, dat$gamma, M)/PBar(0, dat$alpha, dat$beta, dat$gamma, M)
#				#the bh fit found
#				Fs = FMsy(fit$alpha, 2, M)
#				xiHat   = Fs/M 
#				zetaHat = PBar(Fs, fit$alpha, fit$beta, fit$gamma, M)/PBar(0, fit$alpha, fit$beta, fit$gamma, M) #xiHat/(2*xiHat+1) #
#				md = optimize(myDist, c(0, max(xiSims)), xi=dat$xi, zeta=dat$zeta)$objective				
#
#				#NOTE: replace with a propper observation uncertainty
#                                if(length(fit$rsCov)==0){ v=0; lbv=0 }else{ v=getlFV(fit)$lFV; lbv=fit$rsCov['lbeta', 'lbeta'] }
#                                #print( c(dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), v) )
#                                D[di,] = c(dat$xi, dat$zeta, xiInv, zetaInv, xiHat, zetaHat, md, log(Fs), v, fit$lbeta, lbv)
#				#print(dim(D))
#				di = di+1
#			}
#		}
#	}
#	#
#	return(D)
#}
#
#zetaSims = seq(0.15, 0.7, 0.1) #0.01 #seq(0.1, 0.8, 0.05) 	#rev(seq(0.1, 0.80, 0.01)) #
#xiSims =   seq(0, 4.5, 0.5) #0.05 #seq(0.5, 3.5, 0.25)		#seq(0.5, 3.5, 0.05)       #
##zetaSims = seq(0.15, 0.9, 0.025)
##xiSims =   seq(0, 4.5, 0.25) #0.
#
