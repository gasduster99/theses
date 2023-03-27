rm(list=ls())

#
library(VGAM)
library(boot)
library(gMOIP)
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
pellaSRR = function(P, alpha, beta, gamma){ alpha*P/(gamma-1)*(1-(P/beta)^(gamma-1)) }

#
PBar = function(ff, alpha, beta, gamma){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }

#
FMsy = function(alpha, gamma){ alpha/gamma }

#
getPar = function(ff, zeta){
        #
        gamma = c()
        for(z in zeta){
                gammas = c(lambertW0(z*log(z))/log(z), lambertWm1(z*log(z))/log(z))
                gamma = c(gamma, gammas[1+(z>(1/exp(1)))])
        }
        alpha = ff*gamma
        #
        return(cbind(alpha, gamma))
}

#
myDist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-(1/2))^2) }

#
getlFV = function(fit, MM=10^4, samples=F){
        #
        who = c('lalpha')
        C = fit$rsCov[who, who]
        m = c(fit$lalpha)
        sam = rnorm(MM, m, sqrt(C)) #rmvnorm(M, m, C)
        #
        als = exp(sam)
        lFs = sapply(als, function(a){log(FMsy(a, fit$gamma))})
        #
        out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T))
        if( samples ){ out$lFSamples=lFs }
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
                if( dat$xi<xiRange[1] | dat$xi>xiRange[2] ){ next }

                #
                xiBin   = strsplit(f, "_")[[1]]
                zetaBin = xiBin[3]
                xiBin   = xiBin[2]
                #
                zetaBin = gsub("zeta", "", zetaBin)
                zetaBin = as.numeric(gsub(".rda", "", zetaBin))
                xiBin = as.numeric(gsub("xi", "", xiBin))
		
		#
                Fs = FMsy(fit$alpha, fit$gamma)
                xiHat = Fs #/M
                #
                BStar = PBar(Fs, fit$alpha, fit$beta, fit$gamma)
                BZero = PBar(0, fit$alpha, fit$beta, fit$gamma)
                zetaHat = BStar/BZero
		
		#
                md = optimize(myDist, c(0, max(xiRange)), xi=dat$xi, zeta=dat$zeta)$objective

                #
                if(length(fit$rsCov)==0){ lfv=0; lbv=0 }else{
                        #
			lfv=getlFV(fit)$lFV
			lbv=fit$rsCov['lbeta', 'lbeta']
                }
		
		#
                D[i,] = c(xiBin, zetaBin, dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), lfv, log(BZero), lbv)
                i = i+1
        }
        #
        return(D)
}

#
#DATA
#

#
P0 = 10000
mod = "FlatT30" #"ExpT45" #
place = sprintf("./modsPT%s/", mod)

#
xiRes = 0.5
zetaTop = 0.6 #0.7
zetaBot = 0.3#0.2 #0.1
xiBot = 0.1 #0.5
xiTop = 0.7 #3.5

#
f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda"))[1] )

#
D = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
D = D[D$lFV>0 & D$lKV>0,]
D = D[!(round(D$xiHat,1)<0.1 & D$zetaSeed>0.5),] #D[!(D$zetaSeed>0.6),] #D[!(D$xiHat<=D$xiSeed*3/4 & D$zetaSeed>0.5),]#

#
#GP INTERPOLATION
#

#First fit the lF Model and make lF related Predictions

#pick a polynomial mean function
lFy = D$lF
lFV = diag(D$lFV)
lFX = cbind(1, D$xiSeed, D$zetaSeed)

#
xAug = seq(xiBot, xiTop, 0.01) #seq(0.5, 4, 0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
aug = cbind(rep(1, length(xAug)), xAug, 1/(2))
lFX = rbind(lFX, aug)
lFy = c(lFy, log(xAug))
lFV = diag(c(D$lFV, rep(mean(D$lFV)*10^-10, length(xAug)))) #mean(D$lFV)
#lFV = diag( apply(cbind((diag(lFV)), 1e-2), 1, max)*10 )

#
lFaxes = lFX[,2:3]
registerData(lFy, lFX, lFaxes, lFV)
par = c(l1=0.5, l2=0.5, s2=0.5)
lFFit = gpMAP(par, hessian=F, psiSample=F, lower=c(rep(eps(), 3))) #lower=c(0.3, rep(eps(), 2)))
writeLines("lF Fit:")
print(lFFit)

#lF prediction
zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.005) #length.out=)     #seq(0.15, 0.35, 0.001)  #rev(seq(0.1, 0.
xiStar   = seq(min(D$xiSeed), max(D$xiSeed), 0.01) #length.out=)          #seq(1, 3.5, 0.005)
lFXStar = cbind(1, expand.grid(xiStar, zetaStar))
#
vchull = D[chull(D[,c('xiBin', 'zetaBin')]),c('xiBin', 'zetaBin')]
mask = inHull(lFXStar[,c(2,3)], vchull)>=0 
#mask = sapply(1:nrow(lFXStar), function(i){
#                #
#                xi = lFXStar[i,2]
#                zeta = lFXStar[i,3]
#                #
#                bot = min(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#                top = max(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#                #
#                return( zeta>bot & zeta<top & zeta>zetaBot & zeta<zetaTop )
#                #return(T)
#        }
#)
lFPred = gpPredict(lFXStar, lFXStar[,2:3], lFFit)
lFPred[!mask] = NA

#bias
xiHat = exp(lFPred)
xiBias = sweep(xiHat, 1, xiStar)
#zetaBias = sweep(1/(xiHat+2), 2, zetaStar)
zetaBias = sweep(matrix(0.5, nrow(xiHat), ncol(xiHat)), 2, zetaStar)
zetaBias[!mask] = NA

#
eucBias = mcmapply(function(xiHat, xi, zeta){
                myDist(xiHat, xi, zeta)
        }, xiHat, lFXStar[,2], lFXStar[,3], mc.cores=6 #detectCores()
)
eucBias = matrix(eucBias, nrow=length(xiStar), ncol=length(zetaStar))

#Now fit the lK Model and make lK related Predictions

#
lKy = D$lK
lKV = diag(c(D$lKV))#*10)
lKX = cbind(1, D$xiSeed, D$zetaSeed)

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
kBias = exp(lKPred)-P0
#zetaHat*B0Hat
bMSYHat = exp(lKPred)/2 #(xiHat+2)
bMSYBias = sweep(bMSYHat, 2, zetaStar*P0)

#MSY bias
bMSYStar = zetaStar*P0
fMSYStar = xiStar #*M

#
#PLOT
#

#F* bias

#
freq = c(T,F,F,F,F,F)

#
png(sprintf("fMSYBiasPT%s.png", mod))
nCols = 100
maxBias = abs(max(xiBias, na.rm=T))
minBias = abs(min(xiBias, na.rm=T))
posCols = hcl.colors(round(nCols*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
negCols = hcl.colors(round(nCols*minBias/(maxBias+minBias)), "Blues 2", rev=F)
xCols = c(negCols, "#FFFFFF", posCols)
#
xiMask = xiStar>xiBot & xiStar<xiTop
zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
#
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, xiBias,
        col  = adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        cex.lab = 1.5,
        cex.main= 1.5
)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(xCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.1f", rev(seq(min(xiBias[xiMask, zetaMask], na.rm=T)*M, max(xiBias[xiMask, zetaMask], na.rm=T)*M, length.out=le
#        fill = rev(xCols[show]), #colMap[c(1, 10, 20)], 
#        xpd = NA
#)
dev.off()

#F* relative bias

#
chaosThresh = 3.56995
rSurf = exp(lFPred)*2

#
png(sprintf("fMSYRelBiasPT%s.png", mod))
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
        xlab = TeX("$F_{MSY}$"),#'Xi',
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
        xlab = TeX("$F_{MSY}$"), #'Xi',
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
        xlab = TeX("$F_{MSY}$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiBot),
        zlim = c(-maxAbsXBias, maxAbsXBias),
        add  = T
)
image(xiStar, zetaStar, xiBias/xiStar,
        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(1, max(xiBias/xiStar, 1, na.rm=T)),
        add  = T
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(xCols), length.out=nCols) #20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))), 
#        fill = rev(xCols[show]), 
#        xpd = NA
#)
dev.off()

#zeta bias

#
png(sprintf("zetaBiasPT%s.png", mod))
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
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Bias in Estimated $B_{MSY}/B_0$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        cex.lab = 1.5,
        cex.main= 1.5
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(yCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(zetaBias[xiMask, zetaMask], na.rm=T), max(zetaBias[xiMask, zetaMask], na.rm=T), length.out=le
#        fill = rev(yCols[show]), 
#        xpd = NA
#)
points(D$xiSeed, D$zetaSeed)
dev.off()

#euc bias

#
png(sprintf("directionalBiasPT%s.png", mod))
eucCols = hcl.colors(41, "Reds 2", rev=T)
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, eucBias,
        col  = adjustcolor(eucCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
	zlim = c(0, 2),
        cex.lab = 1.5,
        cex.main= 1.5
)
#curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=0.5)
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
#        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
#        fill = rev(eucCols[show]), 
#        xpd = NA
#)
#points(D$xiSeed, D$zetaSeed)
dev.off()

#
eg = expand.grid(xiStar, zetaStar)
ds = matrix(sqrt(eg[,1]^2+eg[,2]^2), nrow=length(xiStar), ncol=length(zetaStar))
#
png(sprintf("directionalBiasPercentPT%s.png", mod))
eucCols = hcl.colors(41, "Reds 2", rev=T)
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, eucBias/ds,
        col  = adjustcolor(eucCols, alpha.f=0.99),
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
	zlim = c(0, 0.6), #1), #1.5),
        cex.lab = 1.5,
        cex.main= 1.5
)
#curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
thin = c(T,rep(F,165))#125))#135))
quiver(
        lFXStar[w,2][thin], lFXStar[w,3][thin],
        xiBias[w][thin], zetaBias[w][thin],
        scale=0.065, length=0.175, lwd=2
)
show = seq(1, length(eucCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
#        fill = rev(eucCols[show]), 
#        xpd = NA
#)
#points(D$xiSeed, D$zetaSeed)
dev.off()

##
#eg = expand.grid(xiStar, zetaStar)
#ds = matrix(sqrt(eg[,1]^2+eg[,2]^2), nrow=length(xiStar), ncol=length(zetaStar))
ms = matrix(abs(0.5-eg[,2]), nrow=length(xiStar), ncol=length(zetaStar))
#
png(sprintf("directionalBiasPercentSmallPT%s.png", mod))
eucCols = hcl.colors(41, "Reds 2", rev=T)
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 4, 4, 4)+0.1)
image(xiStar, zetaStar, eucBias/ms,
        col  = adjustcolor(eucCols, alpha.f=0.99),
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
	#zlim = c(0, 0.6), #1), #1.5),
        zlim = c(1, 3.5), #22.3),
	cex.lab = 1.5,
        cex.main= 1.5
)
#curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
thin = c(T,rep(F,165))#125))#135))
quiver(
        lFXStar[w,2][thin], lFXStar[w,3][thin],
        xiBias[w][thin], zetaBias[w][thin],
        scale=0.065, length=0.175, lwd=2
)
show = seq(1, length(eucCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
#        fill = rev(eucCols[show]), 
#        xpd = NA
#)
#points(D$xiSeed, D$zetaSeed)
dev.off()

#K bias

#
png(sprintf("kBiasPT%s.png", mod))
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
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(kCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"),
#        sprintf("%1.0f", round(rev(seq(min(kBias[xiMask, zetaMask], na.rm=T)*M, max(kBias[xiMask, zetaMask], na.rm=T)*M, length.ou
#        fill = rev(kCols[show]),
#        xpd = NA
#)
dev.off()

#K relative bias

#
png(sprintf("kRelBiasPT%s.png", mod))
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
image(xiStar, zetaStar, kBias/P0,
        col  = adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $B_0$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(-maxAbsXBias, maxAbsXBias),
        cex.lab = 1.5,
        cex.main= 1.5
)
contour(xiStar[15:16], zetaStar, (kBias/P0)[15:16,],
        method = "simple",
        levels = seq(-1, 1, 0.1),
        labcex = 1.025,
        lwd  = eps(),
        xlab = TeX("$F_{MSY}$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiBot),
        zlim = c(-maxAbsXBias, maxAbsXBias),
        add  = T
)
end = rev(length(xiStar)-c(15:16))
contour(xiStar[end], zetaStar, (kBias/P0)[end,],
        method = "simple",
        levels = seq(-1, 1, 0.1),
        labcex = 1.025,
        lwd  = eps(),
        xlab = TeX("$F_{MSY}$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiBot),
        zlim = c(-maxAbsXBias, maxAbsXBias),
        add  = T
)
image(xiStar, zetaStar, kBias/P0,
        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(1, max(kBias/P0, 1, na.rm=T)),
        add  = T
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(xCols), length.out=nCol) #20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
#        fill = rev(xCols[show]), 
#        xpd = NA
#)
dev.off()

#B* bias

#
png(sprintf("bMSYBiasPT%s.png", mod))
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
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Bias in Estimated $B_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        cex.lab = 1.5,
        cex.main= 1.5
)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(bMSYCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"), 
#        sprintf("%1.0f", round(rev(seq(min(bMSYBias[xiMask, zetaMask], na.rm=T)*M, max(bMSYBias[xiMask, zetaMask], na.rm=T)*M, len
#        fill = rev(bMSYCols[show]),
#        xpd = NA
#)
dev.off()

#B* relative bias

#
png(sprintf("bMSYRelBiasPT%s.png", mod))
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
image(xiStar, zetaStar, t(t(bMSYBias)/(P0*zetaStar)),
        col  = adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $B_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop)
        ,zlim = c(-maxAbsXBias, maxAbsXBias),
        cex.lab = 1.5,
        cex.main= 1.5
)
contour(xiStar[15:16], zetaStar, (t(t(bMSYBias)/(P0*zetaStar)))[15:16,],
        method = "simple",
        levels = seq(-1, 1, 0.1),
        labcex = 1.025,
        lwd  = eps(),
        xlab = TeX("$F_{MSY}$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiBot),
        zlim = c(-maxAbsXBias, maxAbsXBias),
        add  = T
)
end = rev(length(xiStar)-c(15:16))
contour(xiStar[end], zetaStar, (t(t(bMSYBias)/(P0*zetaStar)))[end,],
        method = "simple",
        levels = seq(-1, 1, 0.1),
        labcex = 1.025,
        lwd  = eps(),
        xlab = TeX("$F_{MSY}$"), #'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiBot),
        zlim = c(-maxAbsXBias, maxAbsXBias),
        add  = T
)
image(xiStar, zetaStar, t(t(bMSYBias)/(P0*zetaStar)),
        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(1, max(t(t(bMSYBias)/(P0*zetaStar)), 1, na.rm=T)),
        add  = T
)
#points(D$xiSeed, D$zetaSeed)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
show = seq(1, length(xCols), length.out=nCol) #20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
#        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
#        fill = rev(xCols[show]),
#        xpd = NA
#)
dev.off()


