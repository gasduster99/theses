rm(list=ls())

#
library(gMOIP)
library(dismo)
library(pracma)
library(mvtnorm)
library(rootSolve)
library(latex2exp)
library(matrixStats)
library(RColorBrewer)
#
source('prodClass0.1.1.r')
source('gpFunkGaussLxLy.r')

#
#FUNCTIONS
#

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

##PT
#SRR = function(P, lalpha, lbeta, gamma){ 
#	#
#        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
#}
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
#diff = function(f, x, h=0.00001){
#	(f(x)-f(x-h))/h
#}
#
##
#dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
#        #linearly interpolate catches
#        ft = floor(t)
#        q  = (t-ft)
#        Cl = catch[ft]
#        Cu = catch[ft+1]
#        C = q*Cu + (1-q)*Cl
#        if(q==0){ C=Cl }
#        #
#        Fmsy = exp(lalpha)/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
#        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
#        out = R - Fmsy*C*P
#        #
#        return( list(out) )
#}


##
#SRR = function(B, alpha, beta, gamma){
#        alpha*B*(1-beta*gamma*B)^(1/gamma)
#}
#
SRR = function(B, lalpha, lbeta, gamma){
        exp(lalpha)*B*(1-exp(lbeta)*gamma*B)^(1/gamma)
}

#
yield = function(B, lalpha, lbeta, gamma){
	SRR(B, lalpha, lbeta, gamma)-M*B
}
yield = Vectorize(yield, c("B"))

#
FMsy = function(alpha, gamma, M){
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
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) } #((alpha/(M

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
dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        if( gamma==-1 ){ Fmsy = sqrt(exp(lalpha)*M)-M
        } else {
                #
                FUpper = abs(exp(lalpha)-M)
                Fmsy = uniroot.all(function(ff){ a(gamma, ff, M)-exp(lalpha) }, c(0, FUpper))[1]
        }
        #
        R = exp(lalpha)*P*(1-exp(lbeta)*gamma*P)^(1/gamma)
        out = R - (M+Fmsy*C)*P
        #
        return( list(out) )
}

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
getZeta = function(gamma, ff, M){
        #fgfm = ff*gamma/(ff+M)
        #fgfm/(1 + fgfm - (M/ff+M)^gamma)

        #
        (1-((M+ff)/getAlpha(gamma, ff, M))^gamma) / (1-(M/getAlpha(gamma, ff, M))^gamma)
}
z = Vectorize(getZeta, "gamma")

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
                md = stats::optimize(myDist, c(0, max(xiRange)), xi=dat$xi, zeta=dat$zeta)$objective

                #
                if(length(fit$rsCov)==0){ lfv=0; lbv=0 }else{
                        #
                        lV = getlV(fit)
                        #
                        lfv = lV$lFV
                        lbv = lV$lKV
                }

                #
                D[i,] = c(xiBin, zetaBin, dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), lfv, log(BZero), lbv)
                i = i+1
        }
        #
        return(D)
}

#
#MAIN
#

##
#dir = "./modsPellaFineQFixReduxP010000" 
#M = 0.2
##
#dir = "./modsPellaFineQFixRFixP010000" 
#M = 0.2

##
#mod = "FlatNoQ"
#dir = sprintf("./modsPella%s/", mod)
#M = 0.2
#P0 = 10000
#

#
#mod = "ExpT45N150MinCon"; contrast=T; fv=1;
#mod = "ExpT45N150Wide"; contrast=T; #"FlatT30N150WideExpHotStart" #"FlatT30N150Wide"
##mod = "HHardFlatT30N150WWide"; contrast=F; fv=1;
mod = "HHardFlatT30N150WWideExtra"; contrast=T; fv=1;
#mod = "HHardFlatT30N150WWideN28"; contrast=F; fv=1;
#mod = "HHardFlatT30N150WWideN56"; contrast=F; fv=1; #nF=1; mod=sprintf("%snF%s", mod, nF); title="150 Samples"; nar=975;
#mod = "HHardFlatT30N150WWideN84"; contrast=F; fv=1;
#mod = "HHardFlatT30N150WWideN112"; contrast=F;
#directionalBiasSchnuteAnimateSinkHHardFlatT30N150WWideN112X3Z0.35.png
#
dir = sprintf("../gpBias/modsSchnute%s/", mod)
P0 = 10000
M = 0.2

#LHS boundaries
xiMin = 0.5
xiMax = 3.5 #2 #3.5
zetaMin = 0.2
zetaMax = 0.7

##
#datFiles = sprintf("%s%s", dir, list.files(path=dir, pattern=glob2rx("datGen*.rda")))
##
#n = length(datFiles)
#minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
#binTrk  = ceiling(abs(log10(minDiff)))
#
##
##datFiles = sprintf("%s%s", dir, list.files(path=dir, pattern=glob2rx("datGen*.rda")))
#xis = unlist(sapply(datFiles, function(fn){
#		dat = readRDS(fn)
#		if( length(dat$zeta)!=0 & length(dat$xi)!=0){
#			return(dat$xi)
#		}
#	})
#)
#zetas = unlist(sapply(datFiles, function(fn){
#		dat = readRDS(fn)
#		if( length(dat$zeta)!=0 & length(dat$xi)!=0){
#                        return(dat$zeta)
#                }
#	})
#)
#
##
#png(sprintf('%sDesign.png', mod))
#plot(xis, zetas)
#dev.off()
#
##
##TARGET YEILD CURVES
##
#
##directionalBiasSchnuteAnimateSinkHHardFlatT30N150WWideN112X3Z0.35.png
##minimize target norm
#xiTar   = 3    #2.5  #3    ##2.5  #3    #0.5   #3.4
#zetaTar = 0.35 #0.45 #0.55 ##0.45 #0.275 #0.55 #0.275 #0.55
##xiTar   = 3 #1 #3#3.4
##zetaTar = 0.45 #0.275#0.55
#norms = sqrt((xiTar-xis)^2 + (zetaTar-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
##
#fileDat = names(who)
#fileFit = gsub("datGen", "fit", fileDat)
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#datLow = dat
#fitLow = fit
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
##grd = 0:10000
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)
#
##
#png(sprintf("yeildCurveCompare%sX%sZ%s.png", mod, xiTar, zetaTar))
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	xlim=c(0, xMax),
#	ylab="Yield",
#	xlab="B",
#	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#    	cex.lab = 2.5,
#    	cex.main= 2.5
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveTruth%sX%sZ%s.png", mod, xiTar, zetaTar))
#par(mar=c(5, 5, 4, 4)+0.1)
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	xlim=c(0, xMax),
#	ylab="Yield",
#	xlab="B",
#	main="Schnute Truth", #TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#    	cex.lab = 2.5,
#    	cex.main= 2.5
#)
##polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])), 
##	col=adjustcolor(cols[1], alpha.f=0.2),
##        border=F
##)
##curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##rug(dat$N, lwd=3, ticksize=0.05)
###rug(dat$N, line=0.75)
###rug(fit$N, col=cols[1])
##legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveFit%sX%sZ%s.png", mod, xiTar, zetaTar))
#par(mar=c(5, 5, 4, 4)+0.1)
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, n=1000,
#	lwd=1,
#	lty=3, 
#	ylim=c(0, yMax),
#	xlim=c(0, xMax),
#	ylab="Yield",
#	xlab="B",
#	main="BH Fit", #TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#    	cex.lab = 2.5,
#    	cex.main= 2.5
#	#,col=cols[1]
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Truth", "BH MLE"), col=c("black", cols[1]), lwd=c(1,3), lty=c(3,1))
#dev.off()


##directionalBiasSchnuteAnimateSinkHHardFlatT30N150WWideN112X3Z0.35.png
##minimize target norm
#xiTar   = 3    #2.5  #3    ##2.5  #3    #0.5   #3.4
#zetaTar = 0.55 #0.45 #0.55 ##0.45 #0.275 #0.55 #0.275 #0.55
##xiTar   = 3 #1 #3#3.4
##zetaTar = 0.45 #0.275#0.55
#norms = sqrt((xiTar-xis)^2 + (zetaTar-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
##
#fileDat = names(who)
#fileFit = gsub("datGen", "fit", fileDat)
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#datHi = dat
#fitHi = fit


#
#ARROW PLOT
#

#
place = sprintf("../gpBias/modsSchnute%s/", mod)
#
freq = c(T,F,F,F,F,F)

#
xiRes = 0.5
zetaTop = 0.6 #0.7
zetaBot = 0.2 #0.1
xiBot = 0.5
xiTop = 3.5

#
f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda"))[1] )
M = readRDS(f)$M #0.2

#
D = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
D = D[D$lFV>0 & D$lKV>0,]
#Remove for the chull
#0.540   0.270 0.6291777 0.2877782
#0.534   0.347 0.5442120 0.3539416
D = D[D$xiBin!=0.540 & D$zetaBin!=0.270,]
D = D[D$xiBin!=0.534 & D$zetaBin!=0.347,]

#
#GP INTERPOLATION
#

#First fit the lF Model and make lF related Predictions

#
augV = 10^-8
sigV = 2
#pick a polynomial mean function
lFy = D$lF
lFV = diag(D$lFV)
lFX = cbind(1, D$xiSeed, D$zetaSeed)
#diag(lFV)[lFX[,3]<1/exp(1)] = diag(lFV)[lFX[,3]<1/exp(1)]*sigV
#
xAug = seq(0.5, 4, 0.05) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
aug = cbind(rep(1, length(xAug)), xAug, 1/(xAug+2))
lFX = rbind(lFX, aug)
lFy = c(lFy, log(xAug*M))
#lFV = diag(c(D$lFV, rep(mean(D$lFV), length(xAug))))
lFV = diag(c(D$lFV, rep(mean(D$lFV)*augV, length(xAug))))

#
lFaxes = lFX[,2:3]
registerData(lFy, lFX, lFaxes, lFV)
par = c(l1=0.5, l2=0.5, s2=0.5)
lFFit = gpMAP(par, hessian=F, psiSample=F, lower=c(rep(eps(), 3))) #lower=c(0.3, rep(eps(), 2)))
writeLines("lF Fit:")
print(lFFit)

#lF prediction
zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.005) #length.out=)   #seq(0.15, 0.35, 0.001)  #rev(seq(0.1, 0.80, 0.01)) #
xiStar   = seq(min(D$xiSeed), max(D$xiSeed), 0.01) #length.out=)        #seq(1, 3.5, 0.005)
lFXStar = cbind(1, expand.grid(xiStar, zetaStar))
#
vchull = D[chull(D[,c('xiBin', 'zetaBin')]),c('xiBin', 'zetaBin')]
mask = inHull(lFXStar[,c(2,3)], vchull)>=0 
mask1 = sapply(1:nrow(lFXStar), function(i){
                #
                win = 0.15#6 #0.3
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
zetaBias = sweep(1/(xiHat+2), 2, zetaStar)
#zetaBias = sweep(matrix(0.5, nrow(xiHat), ncol(xiHat)), 2, zetaStar)
zetaBias[!mask] = NA

#
eucBias = mcmapply(function(xiHat, xi, zeta){
                myDist(xiHat, xi, zeta)
        }, xiHat, lFXStar[,2], lFXStar[,3], mc.cores=6 #detectCores()
)
eucBias = matrix(eucBias, nrow=length(xiStar), ncol=length(zetaStar))

##Now fit the lK Model and make lK related Predictions
#
##
#lKy = D$lK
#lKV = diag(c(D$lKV))
#lKX = cbind(1, D$xiSeed, D$zetaSeed)
##diag(lKV)[lKX[,3]<(1/exp(1))] = diag(lKV)[lKX[,3]<(1/exp(1))]*sigV
###
##xAug = seq(0.5, 4, 0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
##aug = cbind(rep(1, length(xAug)), xAug, 1/(xAug+2))
#lKX = rbind(lKX, aug)
#lKy = c(lKy, rep(log(P0), length(xAug)))
##lKV = diag(c(D$lKV, rep(mean(D$lKV), length(xAug))))
#lKV = diag(c(D$lKV, rep(mean(D$lKV)*augV, length(xAug))))
#
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
#kBias = exp(lKPred)-P0
##zetaHat*B0Hat
#bMSYHat = exp(lKPred)/(xiHat+2)
#bMSYBias = sweep(bMSYHat, 2, zetaStar*P0)
#
##MSY bias
#bMSYStar = zetaStar*P0
#fMSYStar = xiStar*M

#extra fill
xiFill   = seq(0, min(D$xiSeed), 0.01)
dotStar = expand.grid(xiFill, zetaStar)

##
#png(sprintf("directionalBiasSchnuteAnimateBase%sX%sZ%s.png", mod, xiTar, zetaTar))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias,
#        col  = adjustcolor(eucCols, alpha.f=0.99),
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop), #c(0, xiTop), # 
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.05, length=0.15 #scale=0.025
#)
#show = seq(1, length(eucCols), length.out=20)
#legend(grconvertX(421.25, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(eucBias, na.rm=T), max(eucBias, na.rm=T), length.out=length(show)))),
#        fill = rev(eucCols[show]), 
#        xpd = NA
#)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#title = c(TeX("Low Contrast"), TeX("High Contrast"))
#png(sprintf("directionalBiasSchnuteMinDist%s.png", mod))
##
#eg = expand.grid(xiStar, zetaStar)
#ms = apply(eg, 1, function(r){ stats::optimize(function(x){ sqrt((r[1]-x)^2 + (r[2]-(1/(x+2)))^2) }, c(0,4))$objective })
#ms = matrix(ms, nrow=length(xiStar), ncol=length(zetaStar))
##mdThresh = 
##low:
##fivenum(log(eucBias/ms))
##[1] -8.8990744 -1.0013181 -0.1577205  0.5057364  1.1217219
##high:
##fivenum(log(eucBias/ms))
##[1] 5.717565e-07 1.252115e+00 1.741566e+00 2.115137e+00 7.808551e+00
##
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ms,
#        col  = adjustcolor(eucCols, alpha.f=0.99),  #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = title[contrast+1], #TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop), # c(0, xiTop), #
#        zlim = (c(1, 15)), #1.71466050), #c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=4, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,475)) #125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2 #scale=0.05, length=0.15 
#)
##points(D$xiSeed, D$zetaSeed)
#dev.off()

#
#PLOT
#

#
nx = 1000
ny = 2*1000
r = raster(matrix(1:(nx*ny), nx, ny), xmn=0.5, xmx=3.5, ymn=0.2, ymx=0.6)

#
#U = D[c(T,F),c("xiSeed", "zetaSeed", "xiHat", "zetaHat")] #cbind(D$xiSeed[c(T,F)], D$zetaSeed)
U = D[,c("xiSeed", "zetaSeed", "xiHat", "zetaHat")]
#U = U[U$zetaSeed<0.6,]
#plot(D$xiSeed[c(T,F)], D$zetaSeed[c(T,F)])
#quiver(D$xiSeed[c(T,F)], D$zetaSeed[c(T,F)], D$xiHat[c(T,F)]-D$xiSeed[c(T,F)], D$zetaHat[c(T,F)]-D$zetaSeed[c(T,F)])
#
v = voronoi(U)#, ext=c(0.5, 0))

##
#eg = expand.grid(xiStar, zetaStar)
#ss = matrix(abs(eg[,2]-0.5), nrow=length(xiStar), ncol=length(zetaStar))
ssThresh = 3.11

#
eucBiasObs = sqrt((U$xiHat-U$xiSeed)^2+(U$zetaHat-U$zetaSeed)^2)
ssObs = (U$zetaSeed-0.5)

#eucBiasObs-ssObs)#
rast = rasterize(v, r, eucBiasObs, ext=extent(0.5, 3.5, 0.2, 0.6)) #abs(U$xiHat-U$xiSeed)) # #(eucBiasObs-U$zetaSeed)) #(eucBias-ss) )#U[,2])
extent(rast) = extent(0.5, 3.5, 0.2, 0.6)
#(eucBias-ss)
eucCols = c(hcl.colors(41, "Reds 2", rev=T))#, "black")
eucCols = adjustcolor(eucCols, alpha.f=0.8)

#
png(sprintf("obsDirectionalBiasSub%s.png", mod))
plot(clamp(rast, 0, ssThresh),
        col=eucCols,
        xlim=c(0.5, 3.5),
        ylim=c(0.2, 0.6),
        asp=NA,
        zlim=c(0, ssThresh),
        ylab="Bmsy/B0",
        xlab="Fmsy/M",
        main="Nearest Neighbor Interpolation of Raw Data"
)
#points(U)
quiver(U$xiSeed, U$zetaSeed, U$xiHat-U$xiSeed, U$zetaHat-U$zetaSeed, scale=0.065, length=0.175, lwd=2)
#segments(U$xiSeed, U$zetaSeed, U$xiHat, U$zetaHat)
#abline(h=0.5, lwd=3)
curve(1/(x+2), lwd=3, add=T)
dev.off()


##
#png(sprintf("obsDirectionalBiasSubLine%s.png", mod))
#plot(clamp(rast, 0, ssThresh), col=eucCols, xlim=c(0.1, 0.7), ylim=c(0.25, 0.6), asp=NA, zlim=c(0, 0.25))
##points(U)
##quiver(U$xiSeed, U$zetaSeed, U$xiHat-U$xiSeed, U$zetaHat-U$zetaSeed, scale=0.025)
#segments(U$xiSeed, U$zetaSeed, U$xiHat, U$zetaHat)
##abline(h=0.5, lwd=3)
#curve(1/(x+2), lwd=3)
#dev.off()







#
title = c(TeX("Low Contrast"), TeX("High Contrast"))
png(sprintf("directionalBiasSchnuteSub%s.png", mod))
#
eg = expand.grid(xiStar, zetaStar)
ms = apply(eg, 1, function(r){ stats::optimize(function(x){ sqrt((r[1]-x)^2 + (r[2]-(1/(x+2)))^2) }, c(0,4))$objective })
ms = matrix(ms, nrow=length(xiStar), ncol=length(zetaStar))
#
msThresh = 3.11
#
eucCols = hcl.colors(41, "Reds 2", rev=T)
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 5, 4, 4)+0.1)
image(xiStar, zetaStar, eucBias-ms,
        col  = adjustcolor(eucCols, alpha.f=0.99),  #eucCols, #
        xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = title[contrast+1], #TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop), # c(0, xiTop), #
        zlim = (c(0, msThresh)), #1.71466050), #c(0, 1), #0.95),
	cex.lab = 1.5,
        cex.main= 1.5
)
greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
image(xiStar, zetaStar, (eucBias-ms),
        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(msThresh, max((eucBias-ms), msThresh, na.rm=T)), #log(c(lmThresh, max(eucBias/ms, mThresh, na.rm=T))),
        add  = T
)
points(dotStar[freq,1], dotStar[freq,2], pch='.')
#curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
curve(1/(x+2), from=0, to=4, lwd=4, add=T)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
thin = c(T,rep(F,475)) #125))#135))
quiver(
        lFXStar[w,2][thin], lFXStar[w,3][thin],
        xiBias[w][thin], zetaBias[w][thin],
        scale=0.065, length=0.175, lwd=2 #scale=0.05, length=0.15 
)
#points(D$xiSeed, D$zetaSeed)
dev.off()

##
#title = c(TeX("Low Contrast"), TeX("High Contrast"))
#png(sprintf("directionalBiasSchnuteSubSpace%s.png", mod))
##
#eg = expand.grid(xiStar, zetaStar)
#ms = apply(eg, 1, function(r){ stats::optimize(function(x){ sqrt((r[1]-x)^2 + (r[2]-(1/(x+2)))^2) }, c(0,4))$objective })
#ms = matrix(ms, nrow=length(xiStar), ncol=length(zetaStar))
##
#msThresh = 3.11
##
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias-ms,
#        col  = adjustcolor(eucCols, alpha.f=0.99),  #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = title[contrast+1], #TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(0, xiTop), #c(xiBot, xiTop), #
#        zlim = (c(0, msThresh)), #1.71466050), #c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
#image(xiStar, zetaStar, (eucBias-ms),
#        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(msThresh, max((eucBias-ms), msThresh, na.rm=T)), #log(c(lmThresh, max(eucBias/ms, mThresh, na.rm=T))),
#        add  = T
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=4, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,475)) #125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2 #scale=0.05, length=0.15 
#)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#title = c(TeX("Low Contrast"), TeX("High Contrast"))
#png(sprintf("directionalBiasSchnuteSubTitle%s.png", mod))
##
#eg = expand.grid(xiStar, zetaStar)
#ms = apply(eg, 1, function(r){ stats::optimize(function(x){ sqrt((r[1]-x)^2 + (r[2]-(1/(x+2)))^2) }, c(0,4))$objective })
#ms = matrix(ms, nrow=length(xiStar), ncol=length(zetaStar))
##
#msThresh = 3.11
##
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias-ms,
#        col  = adjustcolor(eucCols, alpha.f=0.99),  #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop), # c(0, xiTop), #
#        zlim = (c(0, msThresh)), #1.71466050), #c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#text = c("", "High Contrast")
#mtext(text[contrast+1], line=3, cex=1.5)
#greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
#image(xiStar, zetaStar, (eucBias-ms),
#        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(msThresh, max((eucBias-ms), msThresh, na.rm=T)), #log(c(lmThresh, max(eucBias/ms, mThresh, na.rm=T))),
#        add  = T
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=4, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,475)) #125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2 #scale=0.05, length=0.15 
#)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#png(sprintf("legendSubSchnute.png"), width = 150, height = 400)
#show = seq(1, length(eucCols), length.out=20)
#plot.new()
#legend("top", #grconvertX(421.25, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", seq(0, msThresh, length.out=20)), #round(rev(seq(min(eucBias/ds, na.rm=T), max(eucBias/ds, na.rm=T), length.out=length(show))/5)*5*100)), #length(show)))),
#        fill = eucCols[seq(1, 41, 2)], #adjustcolor(rev(eucCols[show]), alpha.f=0.99), 
#        xpd = NA
#)
#dev.off()
#
##
#png(sprintf("directionalBiasSchnuteNorm%s.png", mod))
##
#eg = expand.grid(xiStar, zetaStar)
#ds = matrix(sqrt(eg[,1]^2+eg[,2]^2), nrow=length(xiStar), ncol=length(zetaStar))
##ms = apply(eg, 1, function(r){ stats::optimize(function(x){ sqrt((r[1]-x)^2 + (r[2]-(1/(x+2)))^2) }, c(0,4))$objective })
##ms = matrix(ms, nrow=length(xiStar), ncol=length(zetaStar))
##
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ds,
#        col  = adjustcolor(eucCols, alpha.f=0.99),  #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = title[contrast+1], #TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop), # c(0, xiTop), #
#        #zlim = (c(1, 15)), #1.71466050), #c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=4, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,475)) #125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2 #scale=0.05, length=0.15 
#)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#png(sprintf("directionalBiasSchnuteAnimatePrecent%sX%sZ%s.png", mod, xiTar, zetaTar))
##
#eg = expand.grid(xiStar, zetaStar)
#ds = matrix(sqrt(eg[,1]^2+eg[,2]^2), nrow=length(xiStar), ncol=length(zetaStar))
##
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ds,
#        col  = adjustcolor(eucCols, alpha.f=0.99),  #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop), # c(0, xiTop), #
#        zlim = c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.05, length=0.15 #scale=0.025
#)
#show = seq(1, length(eucCols), length.out=20)
#legend(grconvertX(421.25, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%s %%", seq(95, 0, -5)), #round(rev(seq(min(eucBias/ds, na.rm=T), max(eucBias/ds, na.rm=T), length.out=length(show))/5)*5*100)), #length(show)))),
#        fill = adjustcolor(rev(eucCols[show]), alpha.f=0.99), 
#        xpd = NA
#)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#png(sprintf("directionalBiasSchnuteAnimateSource%sX%sZ%s.png", mod, xiTar, zetaTar))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ds, #eucBias,
#        col  = adjustcolor(eucCols, alpha.f=0.99), #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop), #c(0, xiTop), #
#        zlim = c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#points(xi, zeta, pch=19, cex=1.5)
#points(xi, zeta, cex=1.5, col='white')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.05, length=0.15 #scale=0.025
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=length(show
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#png(sprintf("directionalBiasSchnuteAnimateSink%sX%sZ%s.png", mod, xiTar, zetaTar))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ds, #eucBias,
#        col  = adjustcolor(eucCols, alpha.f=0.99), #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(0, xiTop), #c(xiBot, xiTop), #
#        zlim = c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#points(xi, zeta, pch=19, cex=1.5)
#points(xi, zeta, cex=1.5, col='white')
#points(D[D$xiSeed==xi & D$zetaSeed==zeta,'xiHat'], D[D$xiSeed==xi & D$zetaSeed==zeta,'zetaHat'], pch=19, cex=1.5, col='red')
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.05, length=0.15 #scale=0.035, length=0.13  #0.025
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=length(show
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#png(sprintf("directionalBiasSchnuteAnimateBoth%sX%sZ%s.png", mod, xiTar, zetaTar))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ds, #eucBias,
#        col  = adjustcolor(eucCols, alpha.f=0.99), #eucCols, #
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}/M$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(0, xiTop), #c(xiBot, xiTop), #
#        zlim = c(0, 1), #0.95),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(dotStar[freq,1], dotStar[freq,2], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##
#points(xi, zeta, pch=19, cex=1.5)
#points(xi, zeta, cex=1.5, col='white')
#points(D[D$xiSeed==xi & D$zetaSeed==zeta,'xiHat'], D[D$xiSeed==xi & D$zetaSeed==zeta,'zetaHat'], pch=19, cex=1.5, col='red')
##
##minimize target norm
#xiTar   = 3    #2.5  #3    ##2.5  #3    #0.5   #3.4
#zetaTar = 0.55 #0.45 #0.55 ##0.45 #0.275 #0.55 #0.275 #0.55
##xiTar   = 3 #1 #3#3.4
##zetaTar = 0.45 #0.275#0.55
#norms = sqrt((xiTar-xis)^2 + (zetaTar-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
#points(xi, zeta, pch=19, cex=1.5)
#points(xi, zeta, cex=1.5, col='white')
#points(D[D$xiSeed==xi & D$zetaSeed==zeta,'xiHat'], D[D$xiSeed==xi & D$zetaSeed==zeta,'zetaHat'], pch=19, cex=1.5, col='red')
##
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.05, length=0.15 #scale=0.035, length=0.13  #0.025
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=length(show
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##zeta bias
#
#
##nCols = 100
##maxBias = abs(max(zetaBias, na.rm=T))
##minBias = abs(min(zetaBias, na.rm=T))
##posCols = hcl.colors(round(nCols*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
##negCols = hcl.colors(round(nCols*minBias/(maxBias+minBias)), "Blues 2", rev=F)
##yCols = c(negCols, "#FFFFFF", posCols)
##
##
#png(sprintf("zetaRelBiasSchnute%s.png", mod))
#nCols = 21 #50*2
#maxAbsXBias = 1 #max(abs(xiBias/xiStar), na.rm=T)
#posCols = hcl.colors(nCols/2, "Reds 2", rev=T)
#negCols = hcl.colors(nCols/2, "Blues 2", rev=F)
#yCols = c(negCols, "#FFFFFF", posCols)
##
##par(mar=c(5, 4, 4, 5)+0.i)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, t(t(zetaBias)/zetaStar),
#        col  = adjustcolor(yCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated $B_{MSY}/B_0$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	zlim = c(-1, 1),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#contour(xiStar[15:16], zetaStar, (t(t(zetaBias)/zetaStar))[15:16,],
#        method = "simple",
#        levels = seq(-1, 1, 0.1),
#        labcex = 1.025,
#        lwd  = eps(),
#        #xlab = TeX("$F_{MSY}/M$"), #'Xi',
#        #ylab = TeX('$B_{MSY}/B_0$'),
#        #main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiBot),
#        zlim = c(-1,1), #c(-maxAbsXBias, maxAbsXBias),
#        add  = T
#)
#end = rev(length(xiStar)-c(15:16))
#contour(xiStar[end], zetaStar, (t(t(zetaBias)/zetaStar))[end,],
#        method = "simple",
#        levels = seq(-1, 1, 0.1),
#        labcex = 1.025,
#        lwd  = eps(),
#        #xlab = TeX("$F_{MSY}/M$"), #'Xi',
#        #ylab = TeX('$B_{MSY}/B_0$'),
#        #main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiBot),
#        zlim = c(-1, 1), #c(-maxAbsXBias, maxAbsXBias),
#        add  = T
#)
#image(xiStar, zetaStar, t(t(zetaBias)/zetaStar),
#        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
#        #xlab = TeX("$F_{MSY}/M$"),#'Xi',
#        #ylab = TeX('$B_{MSY}/B_0$'),
#        #main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(1, max(zetaBias/zetaStar, 1, na.rm=T)),
#        add  = T
#)
##points(D$xiSeed, D$zetaSeed)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#show = seq(1, length(yCols), length.out=20)
#dev.off()
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(zetaBias[xiMask, zetaMask], na.rm=T), max(zetaBias[xiMask, zetaMask], na.rm=T), length.out=length(show)))),
##        fill = rev(yCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
##dev.off()
#
##
#png(sprintf("fMSYRelBiasSchnute%s.png", mod))
#nCols = 21 #50*2
#maxAbsXBias = 1 #max(abs(xiBias/xiStar), na.rm=T)
#posCols = hcl.colors(nCols/2, "Reds 2", rev=T)
#negCols = hcl.colors(nCols/2, "Blues 2", rev=F)
#xCols = c(negCols, "#FFFFFF", posCols)
##
#xiMask = xiStar>xiBot & xiStar<xiTop
#zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
##
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, xiBias/xiStar,
#        col  = adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}/M$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated $F_{MSY}/M$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(-maxAbsXBias, maxAbsXBias),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#contour(xiStar[15:16], zetaStar, (xiBias/xiStar)[15:16,],
#        method = "simple",
#        levels = seq(-1, 1, 0.1),
#        labcex = 1.025,
#        lwd  = eps(),
#        #xlab = TeX("$F_{MSY}/M$"), #'Xi',
#        #ylab = TeX('$B_{MSY}/B_0$'),
#        #main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiBot),
#        zlim = c(-maxAbsXBias, maxAbsXBias),
#        add  = T
#)
#end = rev(length(xiStar)-c(15:16))
#contour(xiStar[end], zetaStar, (xiBias/xiStar)[end,],
#        method = "simple",
#        levels = seq(-1, 1, 0.1),
#        labcex = 1.025,
#        lwd  = eps(),
#        #xlab = TeX("$F_{MSY}/M$"), #'Xi',
#        #ylab = TeX('$B_{MSY}/B_0$'),
#        #main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiBot),
#        zlim = c(-maxAbsXBias, maxAbsXBias),
#        add  = T
#)
#image(xiStar, zetaStar, xiBias/xiStar,
#        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
#        #xlab = TeX("$F_{MSY}/M$"),#'Xi',
#        #ylab = TeX('$B_{MSY}/B_0$'),
#        #main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(1, max(xiBias/xiStar, 1, na.rm=T)),
#        add  = T
#)
##points(D$xiSeed, D$zetaSeed)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#show = seq(1, length(xCols), length.out=nCols) #20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), 
##        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))), 
##        fill = rev(xCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()



















###
##png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##step=-10
##quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
##	step, step*diffDat,
##	scale=100,
##	col=1:length(dat$N) 
##)
##diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
##	step, step*diffFit,
##	scale=100,
##	col=1:length(dat$N) 
##)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
##dev.off()
##
####
###plot(diffDat[-1], diffFit[-1])
###lines(c(-100, 100), c(-100, 100))
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
#NSam = iSam/fit$q
#depSam = NSam/NSam[,1]
#
##
#png(sprintf("bioCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, fit$N, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(NSam, probs=0.025),
#			rev(colQuantiles(NSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#png(sprintf("depCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0)), #colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(depSam, probs=0.025),
#			rev(colQuantiles(depSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#tmp = prodModel$new(
#	dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)}, #mode
#	time=fit$time, catch=fit$catch, M=M,		#BH        #constants
#	alpha=fit$alpha, beta=fit$beta, gamma=-1,	#parameters
#	lalpha=fit$lalpha, lbeta=fit$lbeta,		#reparameterize
#	lq=log(0.00049), lsdo=log(0.01160256)		#nuisance parameters
#)
##tmp = prodModel$new(
##        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #model
##        time=fit$time, catch=fit$catch, M=M,                          #con
##        alpha=fit$alpha, beta=P0, gamma=2, #parameters
##        lalpha=log(fit$alpha), lbeta=log(P0),       #reparameterize
##        lq=log(0.00049), lsdo=log(0.01160256) #log(0.1160256)
##)

##
##dPost = t(mapply( function(a, b){
#dPost = c()
#for(i in 1:nrow(sam)){
#	#
#	aa = sam[i,1]
#	bb = sam[i,2]
#	#
#	tmp$lalpha = aa
#	tmp$lbeta  = bb
#	tmp$iterate('lsode')
#	
#	#
#	out = tmp$N/tmp$N[1]
#	if( length(out)!=length(tmp$time) ){ next } 
#	dPost = rbind(dPost, out)
#
#	#return(numeric(0)) }
#	#return( out )
##}, sam[,1], sam[,2] ))
#}
#
##
#png(sprintf("depPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, colMeans(dPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(dPost, probs=0.025),
#			rev(colQuantiles(dPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N )
##}, sam[,1], sam[,2] ))
#bPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N
#        if( length(out)!=length(tmp$time) ){ next }
#        bPost = rbind(bPost, out)
#}
#
##
#png(sprintf("bioPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N, colQuantiles(bPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, colMeans(bPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(bPost, probs=0.025),
#			rev(colQuantiles(bPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1e\t& %1.3f\t", round(xi, binTrk), dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$round(xi, binTrk)*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", round(xi, binTrk)*M, round(zeta, binTrk), dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#tl = ff*fit$catch*fit$N
#dtl = dff*dat$catch*dat$N
#ltl = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))

########################################################################################################################
##
##MAX, MIN
##
########################################################################################################################
#
##minimize (max, min) norm
#norms = sqrt((xiMax-xis)^2 + (zetaMin-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
###
##xi = 3.5
##zeta = 0.6
##fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
##fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#fileDat = names(who)                            #sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zet
#fileFit = gsub("datGen", "fit", fileDat)
#
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
##grd = 0:fit$N0
##lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
##qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)
#
##
#png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#    	cex.lab = 1.5,
#    	cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#    	cex.lab = 1.5,
#    	cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
#
###
##png(sprintf("curveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##cols = brewer.pal(9, "Set1")
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
##	col=adjustcolor(cols[1], alpha.f=0.2),
##        border=F
##)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##rug(dat$N, lwd=3, ticksize=0.05)
###rug(dat$N, line=0.75)
###rug(fit$N, col=cols[1])
##legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
##dev.off()
##
###
##png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##step=-10
##quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
##	step, step*diffDat,
##	scale=100,
##	col=1:length(dat$N) 
##)
##diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
##	step, step*diffFit,
##	scale=100,
##	col=1:length(dat$N) 
##)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
##dev.off()
##
####
###bioCol = brewer.pal(9, "YlGnBu")[c(-1, -2)]
####plot(diffDat[-1], diffFit[-1], col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
####lines(c(-100, 100), c(-100, 100))
###plot(SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma), SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma), col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
###lines(c(-10000, 10000), c(-10000, 10000))
####
###bioCol = brewer.pal(9, "YlGnBu")[c(-1, -2)]
####plot(diffDat[-1], diffFit[-1], col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
####lines(c(-100, 100), c(-100, 100))
###srrBar = mapply(function(a, b){SRR(fit$N, a, b, fit$gamma)}, sam[,1], sam[,2])
###qSrr = rowQuantiles(srrBar, probs=c(0.025, 0.5, 0.975))
###plot(SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma), SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma), col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
###for(i in 1:nrow(qSrr)){
###	lines(rep(SRR(dat$N[i], dat$lalpha, dat$lbeta, dat$gamma), 2), c(qSrr[i,1], qSrr[i,3]), lwd=3, col=map2color(dat$N[i], bioCol, limits=c(0, 10000)))
###}
###lines(c(-10000, 10000), c(-10000, 10000))
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
#NSam = iSam/fit$q
#depSam = NSam/NSam[,1]
#
##
#png(sprintf("bioCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, fit$N, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(NSam, probs=0.025),
#			rev(colQuantiles(NSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#png(sprintf("depCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(depSam, probs=0.025),
#			rev(colQuantiles(depSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##dPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N/tmp$N[1] )
##}, sam[,1], sam[,2] ))
#dPost = c()
#for(i in 1:nrow(sam)){
#	#
#	aa = sam[i,1]
#	bb = sam[i,2]
#	#
#	tmp$lalpha = aa
#	tmp$lbeta  = bb
#	tmp$iterate('lsode')
#	
#	#
#	out = tmp$N/tmp$N[1]
#	if( length(out)!=length(tmp$time) ){ next } 
#	dPost = rbind(dPost, out)
#
#	#return(numeric(0)) }
#	#return( out )
##}, sam[,1], sam[,2] ))
#}
#
##
#png(sprintf("depPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, colMeans(dPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(dPost, probs=0.025),
#			rev(colQuantiles(dPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N )
##}, sam[,1], sam[,2] ))
#bPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N
#        if( length(out)!=length(tmp$time) ){ next }
#        bPost = rbind(bPost, out)
#}
#
#
##
#png(sprintf("bioPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N, colQuantiles(bPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, colMeans(bPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(bPost, probs=0.025),
#			rev(colQuantiles(bPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#tr = ff*fit$catch*fit$N
#dtr = dff*dat$catch*dat$N
#ltr = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))
#
###
##xi = 0.5 		#2.3 	#
##zeta = 0.2 	#0.25 	#
#
########################################################################################################################
##
##MIN, MIN
##
########################################################################################################################
#
##minimize (min, min) norm
#norms = sqrt((xiMin-xis)^2 + (zetaMin-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
#
##
#fileDat = names(who) 				#sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#fileFit = gsub("datGen", "fit", fileDat)	#sprintf('%s/fit_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
#
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)
#
##
#png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
#
###
###grd = 0:fit$N0
##lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
##qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
##
###
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)
##
###
##png(sprintf("curveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##cols = brewer.pal(9, "Set1")
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
##	col=adjustcolor(cols[1], alpha.f=0.2),
##        border=F
##)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##rug(dat$N, lwd=3, ticksize=0.05)
###rug(dat$N, line=0.75)
###rug(fit$N, col=cols[1])
##legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
##dev.off()
##
###
##png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##step=-10
##quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
##	step, step*diffDat,
##	scale=100,
##	col=1:length(dat$N) 
##)
##diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
##	step, step*diffFit,
##	scale=100,
##	col=1:length(dat$N) 
##)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
##dev.off()
##
#####
####plot(diffDat, diffFit)
####lines(c(-100, 100), c(-100, 100))
####
###bioCol = brewer.pal(9, "YlGnBu")[c(-1, -2)]
####plot(diffDat[-1], diffFit[-1], col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
####lines(c(-100, 100), c(-100, 100))
###srrBar = mapply(function(a, b){SRR(fit$N, a, b, fit$gamma)}, sam[,1], sam[,2])
###qSrr = rowQuantiles(srrBar, probs=c(0.025, 0.5, 0.975))
###plot(SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma), SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma), col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
###for(i in 1:nrow(qSrr)){
###	lines(rep(SRR(dat$N[i], dat$lalpha, dat$lbeta, dat$gamma), 2), c(qSrr[i,1], qSrr[i,3]), lwd=3, col=map2color(dat$N[i], bioCol, limits=c(0, 10000)))
###}
###lines(c(-10000, 10000), c(-10000, 10000))
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
#NSam = iSam/fit$q
#depSam = NSam/NSam[,1]
##depPSam = 
#
##
#png(sprintf("bioCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, fit$N, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(NSam, probs=0.025),
#			rev(colQuantiles(NSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#png(sprintf("depCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(depSam, probs=0.025),
#			rev(colQuantiles(depSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##dPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N/tmp$N[1] )
##}, sam[,1], sam[,2] ))
#dPost = c()
#for(i in 1:nrow(sam)){
#	#
#	aa = sam[i,1]
#	bb = sam[i,2]
#	#
#	tmp$lalpha = aa
#	tmp$lbeta  = bb
#	tmp$iterate('lsode')
#	
#	#
#	out = tmp$N/tmp$N[1]
#	if( length(out)!=length(tmp$time) ){ next } 
#	dPost = rbind(dPost, out)
#
#	#return(numeric(0)) }
#	#return( out )
##}, sam[,1], sam[,2] ))
#}
#
##
#png(sprintf("depPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(dPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, colMeans(dPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(dPost, probs=0.025),
#			rev(colQuantiles(dPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N )
##}, sam[,1], sam[,2] ))
#bPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N
#        if( length(out)!=length(tmp$time) ){ next }
#        bPost = rbind(bPost, out)
#}
#
##
#png(sprintf("bioPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N, colQuantiles(bPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, colMeans(bPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(bPost, probs=0.025),
#			rev(colQuantiles(bPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
#
##(1.0, 0.2)       &              &               &               &               &               &
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )  #ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#bl = ff*fit$catch*fit$N
#dbl = dff*dat$catch*dat$N
#lbl = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))
#
#
########################################################################################################################
##
##MIN, MAX
##
########################################################################################################################
#
#
###
##xi = 3.5
##zeta = 0.2
#
##minimize (min, max) norm
#norms = sqrt((xiMin-xis)^2 + (zetaMax-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
###
##fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
##fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
##
#fileDat = names(who)                            #sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zet
#fileFit = gsub("datGen", "fit", fileDat)
#
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
##grd = 0:xMax
##lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
##qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)
#
##
#png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)
##
###
##png(sprintf("curveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##cols = brewer.pal(9, "Set1")
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
##	col=adjustcolor(cols[1], alpha.f=0.2),
##        border=F
##)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##rug(dat$N, lwd=3, ticksize=0.05)
###rug(dat$N, line=0.75)
###rug(fit$N, col=cols[1])
##legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
##dev.off()
##
###
##png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##step=-10
##quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
##	step, step*diffDat,
##	scale=100,
##	col=1:length(dat$N) 
##)
##diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
##	step, step*diffFit,
##	scale=100,
##	col=1:length(dat$N) 
##)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
##dev.off()
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
#NSam = iSam/fit$q
#depSam = NSam/NSam[,1]
#
##
#png(sprintf("bioCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, fit$N, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(NSam, probs=0.025),
#			rev(colQuantiles(NSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#png(sprintf("depCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(depSam, probs=0.025),
#			rev(colQuantiles(depSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##dPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N/tmp$N[1] )
##}, sam[,1], sam[,2] ))
#dPost = c()
#for(i in 1:nrow(sam)){
#	#
#	aa = sam[i,1]
#	bb = sam[i,2]
#	#
#	tmp$lalpha = aa
#	tmp$lbeta  = bb
#	tmp$iterate('lsode')
#	
#	#
#	out = tmp$N/tmp$N[1]
#	if( length(out)!=length(tmp$time) ){ next } 
#	dPost = rbind(dPost, out)
#
#	#return(numeric(0)) }
#	#return( out )
##}, sam[,1], sam[,2] ))
#}
#
##
#png(sprintf("depPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, colMeans(dPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(dPost, probs=0.025),
#			rev(colQuantiles(dPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N )
##}, sam[,1], sam[,2] ))
#bPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N
#        if( length(out)!=length(tmp$time) ){ next }
#        bPost = rbind(bPost, out)
#}
#
##
#png(sprintf("bioPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N, colQuantiles(bPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, colMeans(bPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(bPost, probs=0.025),
#			rev(colQuantiles(bPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#br = ff*fit$catch*fit$N
#dbr = dff*dat$catch*dat$N
#lbr = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))

########################################################################################################################
##
##MAX, MID HI
##
########################################################################################################################
#
#
###
##xi = 3.5
##zeta = 0.2
#zetaMidHi = 0.45
# 
##minimize norm
#norms = sqrt((xiMax-xis)^2 + (zetaMidHi-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
###
##fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
##fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
##
#fileDat = names(who)                            #sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zet
#fileFit = gsub("datGen", "fit", fileDat)
#
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
##grd = 0:xMax
##lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
##qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)
#
##
#png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)
##
###
##png(sprintf("curveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##cols = brewer.pal(9, "Set1")
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
##	col=adjustcolor(cols[1], alpha.f=0.2),
##        border=F
##)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##rug(dat$N, lwd=3, ticksize=0.05)
###rug(dat$N, line=0.75)
###rug(fit$N, col=cols[1])
##legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
##dev.off()
##
###
##png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
##diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
##	cex.lab = 1.5,
##    cex.main= 1.5
##)
##step=-10
##quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
##	step, step*diffDat,
##	scale=100,
##	col=1:length(dat$N) 
##)
##diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
##	step, step*diffFit,
##	scale=100,
##	col=1:length(dat$N) 
##)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
##dev.off()
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
#NSam = iSam/fit$q
#depSam = NSam/NSam[,1]
#
##
#png(sprintf("bioCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, fit$N, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(NSam, probs=0.025),
#			rev(colQuantiles(NSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#png(sprintf("depCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(depSam, probs=0.025),
#			rev(colQuantiles(depSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##dPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N/tmp$N[1] )
##}, sam[,1], sam[,2] ))
#dPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#        
#        #
#        out = tmp$N/tmp$N[1]
#        if( length(out)!=length(tmp$time) ){ next }
#        dPost = rbind(dPost, out)
#        
#        #return(numeric(0)) }
#        #return( out )
##}, sam[,1], sam[,2] ))
#}
#
##
#png(sprintf("depPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, colMeans(dPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(dPost, probs=0.025),
#			rev(colQuantiles(dPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N )
##}, sam[,1], sam[,2] ))
#bPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N
#        if( length(out)!=length(tmp$time) ){ next }
#        bPost = rbind(bPost, out)
#}
#
#
##
#png(sprintf("bioPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N, colQuantiles(bPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, colMeans(bPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(bPost, probs=0.025),
#			rev(colQuantiles(bPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#br = ff*fit$catch*fit$N
#dbr = dff*dat$catch*dat$N
#lbr = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))
#
########################################################################################################################
##
##MAX, MID LO
##
########################################################################################################################
#
#
##
#zetaMidLo = 0.35
#
##minimize (min, max) norm
#norms = sqrt((xiMax-xis)^2 + (zetaMidLo-zetas)^2)
#who   = which(min(norms)==norms)
#xi    = xis[who]
#zeta  = zetas[who]
##
#fileDat = names(who)                            #sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zet
#fileFit = gsub("datGen", "fit", fileDat)
#
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
##yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
#
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)
#
##
#png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("yeildCurveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#        lwd=3,
#        ylim=c(0, yMax),
#        ylab="Production",
#        xlab="B",
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
#NSam = iSam/fit$q
#depSam = NSam/NSam[,1]
#
##
#png(sprintf("bioCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, fit$N, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(NSam, probs=0.025),
#			rev(colQuantiles(NSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#png(sprintf("depCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(depSam, probs=0.025),
#			rev(colQuantiles(depSam, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
#dPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N/tmp$N[1]
#        if( length(out)!=length(tmp$time) ){ next }
#        dPost = rbind(dPost, out)
#
#        #return(numeric(0)) }
#        #return( out )
##}, sam[,1], sam[,2] ))
#}
##dPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N/tmp$N[1] )
##}, sam[,1], sam[,2] ))
#
##
#png(sprintf("depPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Depletion",
#        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N/dat$N0, lwd=3)
#lines(fit$time, colMeans(dPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(dPost, probs=0.025),
#			rev(colQuantiles(dPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##	#
##	tmp$lalpha = a
##	tmp$lbeta  = b
##	tmp$iterate('lsode')
##	#
##	return( tmp$N )
##}, sam[,1], sam[,2] ))
#bPost = c()
#for(i in 1:nrow(sam)){
#        #
#        aa = sam[i,1]
#        bb = sam[i,2]
#        #
#        tmp$lalpha = aa
#        tmp$lbeta  = bb
#        tmp$iterate('lsode')
#
#        #
#        out = tmp$N
#        if( length(out)!=length(tmp$time) ){ next }
#        bPost = rbind(bPost, out)
#}
#
##
#png(sprintf("bioPostCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#        xlab="Time",
#        ylab="Biomass",
#        ylim=c(min(0, dat$N), max(fit$N, dat$N, colQuantiles(bPost, probs=0.975))),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#lines(dat$time, dat$N, lwd=3)
#lines(fit$time, colMeans(bPost), lwd=3, col='red')
#polygon( c(fit$time, rev(fit$time)),
#                c(
#			colQuantiles(bPost, probs=0.025),
#			rev(colQuantiles(bPost, probs=0.975))
#		),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#br = ff*fit$catch*fit$N
#dbr = dff*dat$catch*dat$N
#lbr = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))





##
#xi = 0.5
#zeta = 0.5
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
#fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
##
#png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#	cex.lab = 1.5,
#    cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##png(sprintf("catch%s.png", mod))
##fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, 1, 0.5)
##dat = readRDS(fileDat)
##plot(dat$time, dat$catch, xlab="Time", ylab="F(t)/F*", ylim=c(0, 2), main="Catch=F* (F(t)/F*) B(t)", type="l", lwd=3)
##dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#
##
#xi = 2
#zeta = 0.5
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
#fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
##
#png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#	cex.lab = 1.5,
#    cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#
##
#xi = 3.5
#zeta = 0.5
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
#fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
##
#png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#	cex.lab = 1.5,
#    cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()
#
###
##png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, xi, zeta))
##diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
##curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
##	lwd=3, 
##	ylim=c(0, yMax),
##	ylab="Production",
##	xlab="B",
##	main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta))
##)
##step=-10
##quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
##	step, step*diffDat,
##	scale=100,
##	col=1:length(dat$N) 
##)
##diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
##curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
##quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
##	step, step*diffFit,
##	scale=100,
##	col=1:length(dat$N) 
##)
##dev.off()
#
##
#png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#dff = FMsy(dat$alpha, dat$gamma, M)
#bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
##
#png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
#plot(-1, -1,
#        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", xi*M, zeta)),
#        xlab="Time",
#        ylab="CPUE",
#        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
#        xlim=c(min(dat$time), max(dat$time)),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#dat$plotMean(add=T)
#fit$plotMean(add=T, col='red')
#fit$plotBand(col='red')
#abline(h=bMsy*dat$q, lty=2)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) ) 
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$alpha, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalpha']), dat$beta, fit$beta, (fit$lbeta-dat$lbeta)/sqrt(C['lbeta', 'lbeta']), fit$sdo) )
#
#
##
#png(sprintf("fts%s.png", mod))
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, 1, 0.5)
#dat = readRDS(fileDat)
##plot(dat$time, dat$catch, xlab="Time", ylab="F(t)/F*", ylim=c(0, 2), main="Catch=F* (F(t)/F*) B(t)", type="l", lwd=3)
#plot(dat$time, dbl, type='l', col=cols[1], lwd=3, xlab="Time", ylab="Catch", ylim=c(min(bl, tl, tr, br), max(bl, tl, tr, br)), main="Catch")
#lines(dat$time, bl, col=cols[1], lty=2, lwd=3)
#lines(dat$time, dtl, col=cols[2], lwd=3)
#lines(dat$time, tl, col=cols[2], lty=2, lwd=3)
#lines(dat$time, dtr, col=cols[3], lwd=3)
#lines(dat$time, tr, col=cols[3], lty=2, lwd=3)
#lines(dat$time, dbr, col=cols[4], lwd=3)
#lines(dat$time, br, col=cols[4], lty=2, lwd=3)
#legend("topright", legend=c(lbl, ltl, ltr, lbr), col=cols, lwd=3)
#dev.off()








##
#png(sprintf("curveCompareX%sZ%s.png", xi, zeta))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta))
#)
#polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()

##
#ff = FMsy(fit$alpha, fit$gamma, M)
#writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))





