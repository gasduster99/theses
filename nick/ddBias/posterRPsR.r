rm(list=ls())

#
library(VGAM)
library(boot)
library(pracma)
library(Ryacas)
library(mvtnorm)
library(graphics)
library(parallel)
library(latex2exp)
library(rootSolve)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)
#
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

#
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) } #((alpha/(M+

#
getBeta = function(alpha, gamma, M, B0){
        (1-(M/alpha)^gamma)/B0/gamma
}
b = Vectorize(getBeta, c("alpha", "gamma"))

#
getAlpha = function(gamma, ff, M){
        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}
a = Vectorize(getAlpha, "gamma")

#
getZeta = function(gamma, ff, M){
        #fgfm = ff*gamma/(ff+M)
        #fgfm/(1 + fgfm - (M/ff+M)^gamma)

        #
        (1-((M+ff)/getAlpha(gamma, ff, M))^gamma) / (1-(M/getAlpha(gamma, ff, M))^gamma)
}
z = Vectorize(getZeta, "gamma")

#
FMsy = function(alpha, gamma, M){
        #
	if(gamma==-1){ return(sqrt(alpha*M)-M) }
	#
        FUpper = alpha-M
        root = uniroot.all(function(ff){ a(gamma, ff, M)-alpha }, c(0, FUpper))

        #root = uniroot.all(function(ff){ 1 - ff/gamma/(M+ff) - (alpha/(M+ff))^(-1/gamma) }, c(0, 
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        ##
        return(root)
}

#
myDist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-(1/(x+2)))^2) }

#
getlFV = function(fit, MM=10^4, samples=F){
	#
	who = c('lalpha')
	C = fit$rsCov[who, who]
	m = c(fit$lalpha)
	sam = rnorm(MM, m, sqrt(C)) #rmvnorm(M, m, C)
	#
	als = exp(sam)
	lFs = sapply(als, function(a){log(FMsy(a, fit$gamma, M))})
	#
	out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T))
	if( samples ){ out$lFSamples=lFs }
	#
	return(out)
}

#
getlV = function(fit, MM=10^4, samples=F){
	#
	who = c('lalpha', 'lbeta')
	C = fit$rsCov[who, who]
	m = c(fit$lalpha, fit$lbeta)
	sam = rmvnorm(MM, m, C) #rnorm(MM, m, sqrt(C)) #
	#
	abs = exp(sam)
	lFs = sapply(abs[,1], function(a){log(FMsy(a, fit$gamma, M))})
	lKs = mapply(function(p1, p2){log(PBar(p1, p2, fit$gamma, 0, M))}, abs[,1], abs[,2])
	#
	out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T), lK=mean(lKs, na.rm=T), lKV=var(lKs, na.rm=T))
	if( samples ){ out$lFSamples=lFs; out$lKSamples=lKs }
	#
	return(out)
}

#
getData = function(dir, xiRange, zetaRange){
	#
	
	#
	fitFiles = sprintf( "%s%s", dir, list.files(path=dir, pattern=glob2rx("fit*.rda")) )
	#
	i = 1
	D = data.frame(
		xiBin=double(), 
		zetaBin=double(), 
		xiSeed=double(), 
		zetaSeed=double(), 
		xiHat=double(), 
		zetaHat=double(), 
		minDist=double(), 
		lF=double(), 
		lFV=double(), 
		lK=double(), 
		lKV=double(), 
		aSeed=double(),
		bSeed=double(),
		gSeed=double(),
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
		Fs = FMsy(fit$alpha, fit$gamma, M)
		xiHat = Fs/M
		#
		BStar = PBar(fit$alpha, fit$beta, fit$gamma, Fs, M)		
		BZero = PBar(fit$alpha, fit$beta, fit$gamma, 0, M)
		zetaHat = BStar/BZero	
		
		#
		md = optimize(myDist, c(0, max(xiRange)), xi=dat$xi, zeta=dat$zeta)$objective
		
		#
		if(length(fit$rsCov)==0){ lfv=0; lbv=0 }else{ 
			#
			lV = getlV(fit)
			#
			lfv = lV$lFV
			lbv = lV$lKV 
		}

		#
		D[i,] = c(xiBin, zetaBin, dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), lfv, log(BZero), lbv, dat$alpha, dat$beta, dat$gamma)	
		i = i+1
	}
	#
	return(D)
}

#
#DELAY
#

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%28W*
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
Fmsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
	#aMin = M*(M+k)/k/W/(1+M*w/k/W)
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10), tol=.Machine$double.eps)$root
}

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

#beta determines Bzero
getBetaDD = function(B0, M, k, w, W, alpha, gamma){
        f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
        out <- tryCatch({
		uniroot(f, c(0, 10), tol=eps())$root
	}, error=function(err){
                out = NA
        })
}

#
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#Logistic:1; Ricker:0.001; BH:-1
getZetaGamma = function(x, M, k, W, aS, a0, gamma){
        #
        w = vbGrow(aS, k, W, a0)
        #
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBetaDD(B0, M, k, w, W, alpha, gamma)
        if(is.na(beta)){return(NaN)}
	#
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = Fmsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaGamma = Vectorize(getZetaGamma, "x")

##
##DATA
##
#
#mod = "HHardFlatT30N150WWideN112" #"HHardExpT45N150M0.05Wide" #"HHardFlatT30N150M0.1WWideN56" # #"HHardExpT45N150M0.3Wide" # #"HHardExpT45N150M0.1Wide" #"ExpT45N150M0.3Wide" #"HHardFlatT30N150WWideExtra" # #"ExpT45N150Wide" ##"ExpT30L4N150Wide" # #"ExpT45N150" #
#place = sprintf("../gpBias/modsSchnute%s/", mod)
#
##
#xiRes = 0.5
#zetaTop = 0.6 #0.7
#zetaBot = 0.2 #0.1
#xiBot = 0.5
#xiTop = 3.5
#
##
#f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda"))[1] )
#M = readRDS(f)$M #0.2
#
##
#D = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
#D = D[D$lFV>0 & D$lKV>0,]
#
##
#zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.005) #length.out=)  #seq(0.15, 0.35, 0.001)  #rev(seq(0.1
#xiStar   = seq(min(D$xiSeed), max(D$xiSeed), 0.01) #length.out=)       #seq(1, 3.5, 0.005)
#lFXStar = cbind(1, expand.grid(xiStar, zetaStar))
#mask = sapply(1:nrow(lFXStar), function(i){
#               #
#               win = 0.3
#               #
#               xi = lFXStar[i,2]
#               zeta = lFXStar[i,3]
#               #step = win*10/0.1
#               bot = c()
#               top = c()
#               for(x in seq(max(xi-win, xiBot), min(xi+win, xiTop), length.out=30)){
#                       #
#                       isWin = x<=(D$xiBin+win) & x>=(D$xiBin-win)
#                       bot = c(bot, min(D$zetaBin[isWin]) )
#                       top = c(top, max(D$zetaBin[isWin]) )
#               }
#               bot = mean(bot)
#               top = mean(top)
#               ##
#               #bot = min(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#               #top = max(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#               #
#               return( zeta>bot & zeta<top & zeta>zetaBot & zeta<zetaTop )
#               #return(T)
#       }
#)

#
#REGIMES
#

#
P0 = 10000
B0 = P0

#
xiRes = 0.5
zetaTop = 0.6 #0.7
zetaBot = 0.2 #0.1
xiBot = 0.5
xiTop = 3.5

#
cols = brewer.pal(9, 'Set1')

#RP Curves
png("rpCurves.png", width=300, height=325)
zetaLine = function(x){getZetaGamma(x, 0.2, 10, 1, 0.1, -1, -1)}
plot(-1, -1, ylim=c(0.17, 0.5), xlim=c(xiBot, xiTop), xlab=TeX('$F_{MSY}/M$'), ylab=TeX('$B_{MSY}/B_0$'), main="Reference Points")
#curve(zetaLine(x), from=0.01, to=4, lwd=3, col=cols[1], add=T)
#zetaLine = function(x){getZetaGamma(x, 0.2, 10, 1, 0.1, -1, 0.001)}
#curve(zetaLine(x), from=0.01, to=4, lwd=6, add=T)
#curve(zetaLine(x), from=0.01, to=4, lwd=3, col=cols[2], add=T)
#zetaLine = function(x){getZetaGamma(x, 0.2, 10, 1, 0.1, -1, -1)}
#curve(zetaLine(x), from=0.01, to=4, lwd=6, add=T)
#curve(zetaLine(x), from=0.01, to=4, lwd=3, col=cols[3], add=T)

#
n=50	#-1.24
gs = seq(-2.69, 1, length.out=n)
#gs==min(abs(gs--1))
gCol = (gs+1)^-0.3#
bhj = which(!is.na(gCol))[1]
for(j in which(is.na(gCol))){
	gCol[j] = gCol[bhj+(bhj-j)]
}
nums = round((gCol-min(gCol))/max(gCol-min(gCol))*100)
cols = sprintf("grey%s", 100-round((gCol-min(gCol))/max(gCol-min(gCol))*100)) #1:100) #sprintf("grey%s", exp(order(abs(gs--1))))
#
i = 1
for(g in gs){
	#
	zetaLine = function(x){
		getZetaGamma(x, 0.2, 10, 1, 0.1, -1, g)
	}
	to = suppressWarnings(stats::optimize(zetaLine, c(0.6, 3.5))$minimum)
	curve(zetaLine(x), from=0.5, to=to, lwd=6, add=T, col=cols[i])
	i = i+1
}
dev.off()

#reg = c(Inf, 1, 0, -1, -Inf)
#
a = 1.5
b = 0.0001
#
b = 1
M = 0.2

#
png('g3Grey.png', width=300, height=325)
curve(SRR(x*B0, a, b, gs[1])/SRR(B0, a, b, gs[1]), 0, 2/b, lwd=3, ylim=c(0, 2.5), ylab=TeX("$R(B)$ / $R(B_0)$"), xlab=TeX("B / $B_0$"), main=TeX("$R(B)=\\alpha B(1-\\beta\\gamma B)^{1/\\gamma}$"), col='white', n=10000)
mtext(TeX("$\\gamma\\in(-2,1)$"), line=0)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
i=n
for(g in rev(gs)){
        #
        #ff = FMsy(a, g, M)
        B0 = PBar(a, b, g, 0, M)
        #a = A(g, ff, M)
        #b = B(a, g, M, B0)
        w = c(3, 1)[round(g, 1)%%1==0]
        curve(SRR(x*B0, a, b, g)/SRR(B0, a, b, g), 0, 5/b, n=10000, col=cols[i], add=T, lwd=4)
        #abline(0, B0*M, col=cols[sum(g<reg)])
        #print( c(g, sum(g<reg)) )
	i = i-1
}
B0 = PBar(a, b, gs[bhj], 0, M)
curve(SRR(x*B0, a, b, gs[bhj])/SRR(B0, a, b, gs[bhj]), 0, 5/b, n=10000, col=cols[bhj], add=T, lwd=3)
#legend("topright", legend=c(TeX(sprintf("$1\\leq\\gamma$")), TeX(sprintf("$%s > \\gamma\\geq %s$",c(1, 0), c(0, -1))), TeX(sprintf("$-1>\\gamma$"))), col=cols[1:4], lwd=3)
dev.off()


