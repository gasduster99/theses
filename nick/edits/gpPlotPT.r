rm(list=ls())

#
library(VGAM)
library(boot)
library(gMOIP)
library(dismo)
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
yStarOnGrid = function(x, xGrid, yGrid){
        #x      : nx2 matrix with columns 'x^1y^0' & 'x^1y^0'
        #xGrid  : mx2 matrix of predicted locations
        #yGrid  : m-vector of predictions over xGrid
        #

        #
        lons = round(x[,1],1) #'x^1y^0'], 1)
        lats = round(x[,2],1) #'x^0y^1'], 1)
        #
        yOut = matrix(NA, nrow=nrow(x), ncol=1)
        for(i in 1:nrow(x)){
                #where = which(xGrid[,'x^1y^0']==lons[i] & xGrid[,'x^0y^1']==lats[i])
                where = which(xGrid[,1]==lons[i] & xGrid[,2]==lats[i])
		print(where)
		yOut[i] = yGrid[where]
        }
        #
        return( yOut )
}

#
#DATA
#

#
nT = 1
nF = 0

#
P0 = 10000
#mod="FlatT30"; contrast=F; fv=1; nF=1; mod=sprintf("%snF%s", mod, nF); title="150 Samples"; nar=975;
#mod="FlatT30"; contrast=F; fv=1; nF=3; mod=sprintf("%snF%s", mod, nF); title="100 Samples"; nar=975;
#mod="FlatT30"; contrast=F; fv=1; nF=6; mod=sprintf("%snF%s", mod, nF); title="50 Samples"; nar=750;
#mod="ExpT45"; contrast=T; fv=1;
mod="ExpT45MinCon2X"; contrast=T; fv=1; title="150 Samples"; nar=850;
#mod="ExpT45MinCon"; contrast=T; fv=1; nF=0; title="100 Samples"; nar=850;
#mod="ExpT45MinCon"; contrast=T; fv=1; nF=1; mod=sprintf("%snF%s", mod, nF); title="50 Samples"; nar=600; #850;
#mod="ExpT45Sig0.1"; contrast=T; fv=10;
##WORKED
#mod="ExpT45Sig0.1N600"; contrast=T; fv=10;
#mod="ExpT45Sig0.15Half"; contrast=T; fv=1;
#mod="ExpT45Sig0.15"; contrast=T; fv=10;
#mod="ExpT45Sig0.15N600"; contrast=T; fv=100;
#mod="ExpT45Sig0.15N900"; contrast=T; fv=100;
#mod="ExpT45Sig0.3"; contrast=T; fv=10;
#mod="ExpT45Sig0.3N600"; contrast=T; fv=10;
#mod="ExpT45Sig0.3N900"; contrast=T; fv=10;
place = sprintf("../ptNew/modsPT%s/", mod)

#
xiRes = 0.5
zetaTop = 0.6 #0.6 #0.7
zetaBot = 0.2 #0.3# #0.1
xiBot = 0.1 #0.5
xiTop = 0.7 #3.5

#
f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda"))[1] )

#
D = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7)) #c(zetaBot, zetaTop)) #
D = D[D$lFV>0 & D$lKV>0,]
##D = D[!D$xiHat>1,]
D = D[!(round(D$xiHat,1)<0.1 & D$zetaSeed>0.5),] #D[!(D$zetaSeed>0.6),] #D[!(D$xiHat<=D$xiSeed*3/4 & D$zetaSeed>0.5),]#
D = D[c(rep(T, nT), rep(F,nF)),]

#
#GP INTERPOLATION
#

#First fit the lF Model and make lF related Predictions

#pick a polynomial mean function
lFy = D$lF
lFV = diag(D$lFV)*fv
lFX = cbind(1, D$xiSeed, D$zetaSeed)

#
xAug = seq(xiBot, xiTop, 0.1) #seq(xiBot, xiTop, 0.1) #seq(xiBot, xiTop, 0.001) #seq(0.5, 4, 0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
aug = cbind(rep(1, length(xAug)), xAug, 1/(2))
lFX = rbind(lFX, aug)
lFy = c(lFy, log(xAug))
lFV = diag(c(D$lFV, rep(mean(D$lFV)*10^-9, length(xAug)))) #mean(D$lFV)
#lFV = diag( apply(cbind((diag(lFV)), 1e-2), 1, max)*10 )

#
lFaxes = lFX[,2:3]
registerData(lFy, lFX, lFaxes, lFV)
par = c(l1=0.5, l2=0.5, s2=0.5)
lFFit = gpMAP(par, hessian=F, psiSample=F, lower=c(rep(eps(), 3))) #lower=c(0.3, rep(eps(), 2)))
writeLines("lF Fit:")
print(lFFit)

#lF prediction
zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.00125) #length.out=)     #seq(0.15, 0.35, 0.001)  #rev(seq(0.1, 0.
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


##
#lKy = D$lK
#lKV = diag(c(D$lKV))#*10)
#lKX = cbind(1, D$xiSeed, D$zetaSeed)
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
#bMSYHat = exp(lKPred)/2 #(xiHat+2)
#bMSYBias = sweep(bMSYHat, 2, zetaStar*P0)
#
##MSY bias
#bMSYStar = zetaStar*P0
#fMSYStar = xiStar #*M


#
#PLOT
#

#
nx = 1000
ny = 2*1000
r = raster(matrix(1:(nx*ny), nx, ny))

#
#U = D[c(T,F),c("xiSeed", "zetaSeed", "xiHat", "zetaHat")] #cbind(D$xiSeed[c(T,F)], D$zetaSeed)
U = D[,c("xiSeed", "zetaSeed", "xiHat", "zetaHat")]
U = U[U$zetaSeed<0.6,]
#plot(D$xiSeed[c(T,F)], D$zetaSeed[c(T,F)])
#quiver(D$xiSeed[c(T,F)], D$zetaSeed[c(T,F)], D$xiHat[c(T,F)]-D$xiSeed[c(T,F)], D$zetaHat[c(T,F)]-D$zetaSeed[c(T,F)])
#
v = voronoi(U)

##
#eg = expand.grid(xiStar, zetaStar)
#ss = matrix(abs(eg[,2]-0.5), nrow=length(xiStar), ncol=length(zetaStar))
ssThresh = 0.25

#
eucBiasObs = sqrt((U$xiHat-U$xiSeed)^2+(U$zetaHat-U$zetaSeed)^2)
ssObs = (U$zetaSeed-0.5)

#eucBiasObs-ssObs)#
rast = rasterize(v, r, abs(U$xiHat-U$xiSeed)) # #(eucBiasObs-U$zetaSeed)) #(eucBias-ss) )#U[,2])
#(eucBias-ss)
eucCols = c(hcl.colors(41, "Reds 2", rev=T))#, "black")
eucCols = adjustcolor(eucCols, alpha.f=0.8)

#
png(sprintf("obsDirectionalBiasSub%s.png", mod))
plot(clamp(rast, 0, ssThresh), 
	col=eucCols, 
	xlim=c(0.1, 0.7), 
	ylim=c(0.2, 0.6), 
	asp=NA, 
	zlim=c(0, ssThresh),
	ylab="Bmsy/B0",
	xlab="Fmsy",
	main="Nearest Neighbor Interpolation of Raw Data"
)
#points(U)
quiver(U$xiSeed, U$zetaSeed, U$xiHat-U$xiSeed, U$zetaHat-U$zetaSeed, scale=0.065, length=0.175, lwd=2)
#segments(U$xiSeed, U$zetaSeed, U$xiHat, U$zetaHat)
abline(h=0.5, lwd=3)
dev.off()

#
png(sprintf("obsDirectionalBiasSubLine%s.png", mod))
plot(clamp(rast, 0, ssThresh), col=eucCols, xlim=c(0.1, 0.7), ylim=c(0.25, 0.6), asp=NA, zlim=c(0, 0.25))
#points(U)
#quiver(U$xiSeed, U$zetaSeed, U$xiHat-U$xiSeed, U$zetaHat-U$zetaSeed, scale=0.025)
segments(U$xiSeed, U$zetaSeed, U$xiHat, U$zetaHat)
abline(h=0.5, lwd=3)
dev.off()

##with color
#
##eucBiasObs-ssObs)#
#rast = rasterize(v, r, U$xiHat-U$xiSeed) # #(eucBiasObs-U$zetaSeed)) #(eucBias-ss) )#U[,2])
##(eucBias-ss)
#redCols = c(hcl.colors(41, "Reds 2", rev=T))#, "black")
#blueCols = c(hcl.colors(41, "Blue 2", rev=T))
#
##
#png(sprintf("obsDirectionalBiasSubLineColor%s.png", mod))
#plot(rast, col=redCols, xlim=c(0.1, 0.7), ylim=c(0.25, 0.6), asp=NA, zlim=c(0, 0.25))
#plot(rast, col=blueCols, xlim=c(0.1, 0.7), ylim=c(0.25, 0.6), asp=NA, zlim=c(-0.25, 0), add=T)
##points(U)
##quiver(U$xiSeed, U$zetaSeed, U$xiHat-U$xiSeed, U$zetaHat-U$zetaSeed, scale=0.025)
#segments(U$xiSeed, U$zetaSeed, U$xiHat, U$zetaHat)
#abline(h=0.5, lwd=3)
#dev.off()

##
#eucBiasObs = mcmapply(function(xiHat, xi, zeta){
#                myDist(xiHat, xi, zeta)
#        }, U$xiHat, U$xiSeed, U$zetaSeed, mc.cores=6 #detectCores()
#)
#eucBiasObs = matrix(eucBiasObs, nrow=length(U$xiSeed), ncol=length(U$zetaSeed))


#F* bias

#
freq = c(T,F,F,F,F,F)

##
#png(sprintf("fMSYBiasPT%s.png", mod))
#nCols = 100
#maxBias = abs(max(xiBias, na.rm=T))
#minBias = abs(min(xiBias, na.rm=T))
#posCols = hcl.colors(round(nCols*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
#negCols = hcl.colors(round(nCols*minBias/(maxBias+minBias)), "Blues 2", rev=F)
#xCols = c(negCols, "#FFFFFF", posCols)
##
#xiMask = xiStar>xiBot & xiStar<xiTop
#zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
##
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, xiBias,
#        col  = adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#show = seq(1, length(xCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.1f", rev(seq(min(xiBias[xiMask, zetaMask], na.rm=T)*M, max(xiBias[xiMask, zetaMask], na.rm=T)*M, length.out=le
##        fill = rev(xCols[show]), #colMap[c(1, 10, 20)], 
##        xpd = NA
##)
#dev.off()
#
##F* relative bias
#
##
#png(sprintf("fMSYRelBiasPT%s.png", mod))
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
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
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
#        xlab = TeX("$F_{MSY}$"), #'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
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
#        xlab = TeX("$F_{MSY}$"), #'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiBot),
#        zlim = c(-maxAbsXBias, maxAbsXBias),
#        add  = T
#)
#image(xiStar, zetaStar, xiBias/xiStar,
#        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(1, max(xiBias/xiStar, 1, na.rm=T)),
#        add  = T
#)
##points(D$xiSeed, D$zetaSeed)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#show = seq(1, length(xCols), length.out=nCols) #20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), 
##        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))), 
##        fill = rev(xCols[show]), 
##        xpd = NA
##)
#dev.off()
#
##zeta bias
#
##
#png(sprintf("zetaBiasPT%s.png", mod))
#nCols = 100
#maxBias = abs(max(zetaBias, na.rm=T))
#minBias = abs(min(zetaBias, na.rm=T))
#posCols = hcl.colors(round(nCols*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
#negCols = hcl.colors(round(nCols*minBias/(maxBias+minBias)), "Blues 2", rev=F)
#yCols = c(negCols, "#FFFFFF", posCols)
##
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, zetaBias,
#        col  = adjustcolor(yCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main = TeX("Bias in Estimated $B_{MSY}/B_0$"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
##points(D$xiSeed, D$zetaSeed)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#show = seq(1, length(yCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(zetaBias[xiMask, zetaMask], na.rm=T), max(zetaBias[xiMask, zetaMask], na.rm=T), length.out=le
##        fill = rev(yCols[show]), 
##        xpd = NA
##)
#points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##euc bias
#
##
#png(sprintf("directionalBiasPT%s.png", mod))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias,
#        col  = adjustcolor(eucCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	zlim = c(0, 2),
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=0.5)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.025
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#

##
#title = c(TeX("Low Contrast"), TeX("High Contrast"))

#
eg = expand.grid(xiStar, zetaStar)
ss = matrix(abs(eg[,2]-0.5), nrow=length(xiStar), ncol=length(zetaStar))
ssThresh = 0.25
#
png(sprintf("directionalBiasSubPT%s.png", mod))
eucCols = hcl.colors(41, "Reds 2", rev=T)
#par(mar=c(5, 4, 4, 5)+0.1)
par(mar=c(5, 5, 4, 4)+0.1)
image(xiStar, zetaStar, (eucBias-ss),
        col  = adjustcolor(eucCols, alpha.f=0.99),
        xlab = TeX("$F_{MSY}$"),
        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
        main = title, #title[contrast+1], #TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
	zlim = c(0, ssThresh), #0.6), #1), #
        cex.lab = 1.5,
        cex.main= 1.5
)
greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
image(xiStar, zetaStar, (eucBias-ss),
        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
        xlab = TeX("$F_{MSY}$"),#'Xi',
        ylab = TeX('$B_{MSY}/B_0$'),
        ylim = c(zetaBot, zetaTop),
        xlim = c(xiBot, xiTop),
        zlim = c(ssThresh, max((eucBias-ss), ssThresh, na.rm=T)), #log(c(lmThresh, max(eucBias/ms, mThresh, na.rm=T))),
        add  = T
)
#curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
abline(h=0.5, lwd=3)
points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
thin = c(T,rep(F,nar)) #975)) #850))#125))#135))
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

png('subLegnd.png', width = 150, height = 400)
plot.new()
legend("top", #grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
        sprintf("%1.2f", seq(0, ssThresh, length.out=20)), #rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), lengt
        fill = eucCols[seq(1, 41, 2)], #rev(eucCols[show]), 
        xpd = NA
)
dev.off()

##
#eg = expand.grid(xiStar, zetaStar)
#ss = matrix(abs(eg[,2]-0.5), nrow=length(xiStar), ncol=length(zetaStar))
##
#png(sprintf("directionalBiasLogSubPT%s.png", mod))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, log(eucBias-ss),
#        col  = adjustcolor(eucCols, alpha.f=0.99),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = title[contrast+1], #TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	zlim = c(-7, 0.42), #0.6), #1), #
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,850))#125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#eg = expand.grid(xiStar, zetaStar)
#ds = matrix(sqrt(eg[,1]^2+eg[,2]^2), nrow=length(xiStar), ncol=length(zetaStar))
##
#png(sprintf("directionalBiasPercentPT%s.png", mod))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ds,
#        col  = adjustcolor(eucCols, alpha.f=0.99),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	zlim = c(0, 1.5), #0.6), #1), #
#        cex.lab = 1.5,
#        cex.main= 1.5
#)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,165))#125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
###
##eg = expand.grid(xiStar, zetaStar)
##ds = matrix(sqrt(eg[,1]^2+eg[,2]^2), nrow=length(xiStar), ncol=length(zetaStar))
#ms = matrix(abs(0.5-eg[,2]), nrow=length(xiStar), ncol=length(zetaStar))
#mThresh = 3.6
##
#png(sprintf("directionalBiasPercentSmallPT%s.png", mod))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, eucBias/ms,
#        col  = adjustcolor(eucCols, alpha.f=0.99),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	#zlim = c(0, 0.6), #1), #1.5),
#        zlim = c(1, mThresh), #22.3),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
#image(xiStar, zetaStar, eucBias/ms,
#        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(mThresh, max(eucBias/ms, mThresh, na.rm=T)),
#        add  = T
#)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,165))#125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
#lmThresh = 1.6 #2 #2.41
##
#png(sprintf("directionalBiasPercentSmallLogPT%s.png", mod))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, log(eucBias/ms),
#        col  = adjustcolor(eucCols, alpha.f=0.99),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	#zlim = c(0, 0.6), #1), #1.5),
#        zlim = c(0, lmThresh), #log(c(1, lmThresh)), #22.3),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
#greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
#image(xiStar, zetaStar, log(eucBias/ms),
#        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
#        xlab = TeX("$F_{MSY}$"),#'Xi',
#        ylab = TeX('$B_{MSY}/B_0$'),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#        zlim = c(lmThresh, max(log(eucBias/ms), lmThresh, na.rm=T)), #log(c(lmThresh, max(eucBias/ms, mThresh, na.rm=T))),
#        add  = T
#)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=3)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,165))#125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
##
##title = c(TeX("Low Contrast\\nBias in ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"), TeX("High Contrast\\nBias in ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"))
#title = c(TeX("Low Contrast"), TeX("High Contrast"))
#
##
#chaosThresh = 3.56995
#rSurf = exp(lFPred)*2
#cThresh = 2.4 #2.371653 #max( log(eucBias/ms)[rSurf<chaosThresh], na.rm=T )
##
#png(sprintf("directionalBiasPercentSmallLogChaosPT%s.png", mod))
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##par(mar=c(5, 4, 4, 5)+0.1)
#par(mar=c(5, 5, 4, 4)+0.1)
#image(xiStar, zetaStar, log(eucBias/ms),
#        col  = adjustcolor(eucCols, alpha.f=0.99),
#        xlab = TeX("$F_{MSY}$"),
#        ylab = TeX('$B_{MSY}/B_0$'), #'Zeta',
#        main = title[contrast+1], ##TeX("Bias Direction for ($F_{MSY}$, $B_{MSY}/B_0$) Jointly"),
#        ylim = c(zetaBot, zetaTop),
#        xlim = c(xiBot, xiTop),
#	#zlim = c(0, 0.6), #1), #1.5),
#        zlim = c(0, cThresh), #log(c(1, lmThresh)), #22.3),
#	cex.lab = 1.5,
#        cex.main= 1.5
#)
##greyRed =  colorRampPalette(c(eucCols[length(eucCols)],'grey10'))(4)[2]
##image(xiStar, zetaStar[zetaStar>0.5], log(eucBias/ms)[,zetaStar>0.5], 
##        col  = greyRed, #eucCols[length(eucCols)], #"grey10", #adjustcolor(xCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),#'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        zlim = c(cThresh, max(log(eucBias/ms), cThresh, na.rm=T)), #log(c(lmThresh, max(eucBias/ms, mThresh, na.rm=T))),
##        add  = T
##)
##curve(x/(2*x+1), from=0, to=12, lwd=3, add=T) 
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#abline(h=0.5, lwd=4)
#points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
#w = T #!mask #& xBias<16 #(XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,850)) #165))#125))#135))
#quiver(
#        lFXStar[w,2][thin], lFXStar[w,3][thin],
#        xiBias[w][thin], zetaBias[w][thin],
#        scale=0.065, length=0.175, lwd=2
#)
#show = seq(1, length(eucCols), length.out=20)
##legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
##        sprintf("%1.2f", rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
##        fill = rev(eucCols[show]), 
##        xpd = NA
##)
##points(D$xiSeed, D$zetaSeed)
#dev.off()
#
#png('minDistLegnd.png', width = 150, height = 400)
#plot.new()
#legend("top", #grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", seq(0, cThresh, length.out=20)), #rev(seq(min(eucBias[xiMask, zetaMask], na.rm=T), max(eucBias[xiMask, zetaMask], na.rm=T), length.out=leng
#        fill = eucCols[seq(1, 41, 2)], #rev(eucCols[show]), 
#        xpd = NA
#)
#dev.off()
#
###K bias
##
###
##png(sprintf("kBiasPT%s.png", mod))
###
##maxBias = abs(max(kBias, na.rm=T))
##minBias = abs(min(kBias, na.rm=T))
##nCol = 25
##posCols = hcl.colors(round(nCol*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
##negCols = hcl.colors(round(nCol*minBias/(maxBias+minBias)), "Blues 2", rev=F)
##kCols = c(negCols, "#FFFFFF", posCols)
###
##xiMask = xiStar>xiBot & xiStar<xiTop
##zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
###
###par(mar=c(5, 4, 4, 5)+0.1)
##par(mar=c(5, 5, 4, 4)+0.1)
##image(xiStar, zetaStar, kBias, #xiStar[xiMask], zetaStar[zetaMask], xBias[xiMask, zetaMask],i
##        col  = adjustcolor(kCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),#'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Bias in Estimated $B_0$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
###curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##abline(h=0.5, lwd=3)
##show = seq(1, length(kCols), length.out=20)
###legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"),
###        sprintf("%1.0f", round(rev(seq(min(kBias[xiMask, zetaMask], na.rm=T)*M, max(kBias[xiMask, zetaMask], na.rm=T)*M, length.ou
###        fill = rev(kCols[show]),
###        xpd = NA
###)
##dev.off()
##
###K relative bias
##
###
##png(sprintf("kRelBiasPT%s.png", mod))
###
##nCol = 21 #15*2
##maxAbsXBias = 1 #max(abs(kBias/P0), na.rm=T)
##posCols = hcl.colors(nCol/2, "Reds 2", rev=T)
##negCols = hcl.colors(nCol/2, "Blues 2", rev=F)
##xCols = c(negCols, "#FFFFFF", posCols)
###
##xiMask = xiStar>xiBot & xiStar<xiTop
##zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
###
###par(mar=c(5, 4, 4, 5)+0.1)
##par(mar=c(5, 5, 4, 4)+0.1)
##image(xiStar, zetaStar, kBias/P0,
##        col  = adjustcolor(xCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $B_0$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        zlim = c(-maxAbsXBias, maxAbsXBias),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##contour(xiStar[15:16], zetaStar, (kBias/P0)[15:16,],
##        method = "simple",
##        levels = seq(-1, 1, 0.1),
##        labcex = 1.025,
##        lwd  = eps(),
##        xlab = TeX("$F_{MSY}$"), #'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiBot),
##        zlim = c(-maxAbsXBias, maxAbsXBias),
##        add  = T
##)
##end = rev(length(xiStar)-c(15:16))
##contour(xiStar[end], zetaStar, (kBias/P0)[end,],
##        method = "simple",
##        levels = seq(-1, 1, 0.1),
##        labcex = 1.025,
##        lwd  = eps(),
##        xlab = TeX("$F_{MSY}$"), #'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiBot),
##        zlim = c(-maxAbsXBias, maxAbsXBias),
##        add  = T
##)
##image(xiStar, zetaStar, kBias/P0,
##        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),#'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        zlim = c(1, max(kBias/P0, 1, na.rm=T)),
##        add  = T
##)
###points(D$xiSeed, D$zetaSeed)
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
###curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##abline(h=0.5, lwd=3)
##show = seq(1, length(xCols), length.out=nCol) #20)
###legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
###        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
###        fill = rev(xCols[show]), 
###        xpd = NA
###)
##dev.off()
##
###B* bias
##
###
##png(sprintf("bMSYBiasPT%s.png", mod))
###
##maxBias = abs(max(bMSYBias, na.rm=T))
##minBias = abs(min(bMSYBias, na.rm=T))
##nCol = 25
##posCols = hcl.colors(round(nCol*maxBias/(maxBias+minBias)), "Reds 2", rev=T)
##negCols = hcl.colors(round(nCol*minBias/(maxBias+minBias)), "Blues 2", rev=F)
##bMSYCols = c(negCols, "#FFFFFF", posCols)
###
##xiMask = xiStar>xiBot & xiStar<xiTop
##zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
###
###par(mar=c(5, 4, 4, 5)+0.1)
##par(mar=c(5, 5, 4, 4)+0.1)
##image(xiStar, zetaStar, bMSYBias,
##        col  = adjustcolor(bMSYCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Bias in Estimated $B_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
###curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##abline(h=0.5, lwd=3)
##show = seq(1, length(bMSYCols), length.out=20)
###legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(150, "device"), 
###        sprintf("%1.0f", round(rev(seq(min(bMSYBias[xiMask, zetaMask], na.rm=T)*M, max(bMSYBias[xiMask, zetaMask], na.rm=T)*M, len
###        fill = rev(bMSYCols[show]),
###        xpd = NA
###)
##dev.off()
##
###B* relative bias
##
###
##png(sprintf("bMSYRelBiasPT%s.png", mod))
###
##nCols = 21 #2*15
##maxAbsXBias = 1 #max(abs(t(bMSYBias)/(P0*zetaStar)), na.rm=T)
##posCols = hcl.colors(nCols/2, "Reds 2", rev=T)
##negCols = hcl.colors(nCols/2, "Blues 2", rev=F)
##xCols = c(negCols, "#FFFFFF", posCols)
###
##xiMask = xiStar>xiBot & xiStar<xiTop
##zetaMask = zetaStar>zetaBot & zetaStar<zetaTop
###
###par(mar=c(5, 4, 4, 5)+0.1)
##par(mar=c(5, 5, 4, 4)+0.1)
##image(xiStar, zetaStar, t(t(bMSYBias)/(P0*zetaStar)),
##        col  = adjustcolor(xCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),#'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $B_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop)
##        ,zlim = c(-maxAbsXBias, maxAbsXBias),
##        cex.lab = 1.5,
##        cex.main= 1.5
##)
##contour(xiStar[15:16], zetaStar, (t(t(bMSYBias)/(P0*zetaStar)))[15:16,],
##        method = "simple",
##        levels = seq(-1, 1, 0.1),
##        labcex = 1.025,
##        lwd  = eps(),
##        xlab = TeX("$F_{MSY}$"), #'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiBot),
##        zlim = c(-maxAbsXBias, maxAbsXBias),
##        add  = T
##)
##end = rev(length(xiStar)-c(15:16))
##contour(xiStar[end], zetaStar, (t(t(bMSYBias)/(P0*zetaStar)))[end,],
##        method = "simple",
##        levels = seq(-1, 1, 0.1),
##        labcex = 1.025,
##        lwd  = eps(),
##        xlab = TeX("$F_{MSY}$"), #'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiBot),
##        zlim = c(-maxAbsXBias, maxAbsXBias),
##        add  = T
##)
##image(xiStar, zetaStar, t(t(bMSYBias)/(P0*zetaStar)),
##        col  = "grey10", #adjustcolor(xCols, alpha.f=0.6),
##        xlab = TeX("$F_{MSY}$"),#'Xi',
##        ylab = TeX('$B_{MSY}/B_0$'),
##        main = TeX("Relative Bias in Estimated $F_{MSY}$"),
##        ylim = c(zetaBot, zetaTop),
##        xlim = c(xiBot, xiTop),
##        zlim = c(1, max(t(t(bMSYBias)/(P0*zetaStar)), 1, na.rm=T)),
##        add  = T
##)
###points(D$xiSeed, D$zetaSeed)
##points(lFXStar[!mask,2][freq], lFXStar[!mask,3][freq], pch='.')
###curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##abline(h=0.5, lwd=3)
##show = seq(1, length(xCols), length.out=nCol) #20)
###legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(580, "device"), grconvertY(120, "device"), 
###        sprintf("%1.2f", rev(seq(-maxAbsXBias, maxAbsXBias, length.out=length(show)))),
###        fill = rev(xCols[show]),
###        xpd = NA
###)
##dev.off()
##
###
###RESIDUALS
###
##
####residuals
###take = X[,3]!=0.5
###yRes = y[take]
###XRes = X[take,]
###gpPredRes = gpPredict(XRes, XRes[,2:3], gpFit, asMat=F)
###png(sprintf('pellaResiduals%s.png', mod))
###plot(XRes[,2], XRes[,3], col=map2color(yRes-gpPredRes, hcl.colors(10, "Blue-Red 3", rev=F)), pch=19)
###dev.off()
###
####
###png(sprintf("pellaRes%sHist2.png", mod))
###hist(yRes-gpPredRes)
###dev.off()
##
####
###pdf('residPlot.pdf')
###gD = as.geodata(cbind(grid[,'x^1y^0'], grid[,'x^0y^1'], res))
###plot(gD, lowess=T)
###dev.off()
##
#####
####pdf("residPlot%s.png", mod)
####gD = as.geodata(cbind(grid[,'x^1y^0'], grid[,'x^0y^1'], res))
###lFPredVec = gpPredict(lFXStar, lFXStar[,2:3], lFFit, asMat=F)
###lFPredOnGrid = yStarOnGrid(lFX[,2:3], lFXStar[,2:3], lFPredVec)
##
