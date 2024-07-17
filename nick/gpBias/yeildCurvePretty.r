rm(list=ls())

#
library(gMOIP)
library(pracma)
library(mvtnorm)
library(plotrix)
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
chop = function(x){ max(x, 0) }
chop = Vectorize(chop)

#
#HEADER
#

#
mod = "HHardFlatT30N150WWideN112"; contrast=F;

#
dir = sprintf("./modsSchnute%s/", mod)
P0 = 10000
M = 0.2

#LHS boundaries
xiMin = 0.5
xiMax = 3.5 #2 #3.5
zetaMin = 0.2
zetaMax = 0.7
#
datFiles = sprintf("%s%s", dir, list.files(path=dir, pattern=glob2rx("datGen*.rda")))
#
n = length(datFiles)
minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
binTrk  = ceiling(abs(log10(minDiff)))

#
#datFiles = sprintf("%s%s", dir, list.files(path=dir, pattern=glob2rx("datGen*.rda")))
xis = unlist(sapply(datFiles, function(fn){
                dat = readRDS(fn)
                if( length(dat$zeta)!=0 & length(dat$xi)!=0){
                        return(dat$xi)
                }
        })
)
zetas = unlist(sapply(datFiles, function(fn){
                dat = readRDS(fn)
                if( length(dat$zeta)!=0 & length(dat$xi)!=0){
                        return(dat$zeta)
                }
        })
)



#######################################################################################################################
#
#MAX, MID HI
#
#######################################################################################################################


##
#xi = 3.5
#zeta = 0.2
zetaMidHi = 0.45
 
#minimize norm
norms = sqrt((xiMax-xis)^2 + (zetaMidHi-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]
##
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#
fileDat = names(who)                            #sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zet
fileFit = gsub("datGen", "fit", fileDat)

#
dat = readRDS(fileDat)
fit = readRDS(fileFit)
fitS = readRDS('fitSingle_xi3.48_zeta0.48.rda')

#
xMax = max(dat$N0, fit$N0, fitS$N0)
grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
CS = fitS$rsCov[who, who]
mS = c(fitS$lalpha, fitS$lbeta)
samS = rmvnorm(MM, mS, CS)

#
#grd = 0:xMax
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)

#
png(sprintf("yeildCurveCompare%sPrettyX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
x = seq(0, xMax, length.out=1000)
y = yield(x, dat$lalpha, dat$lbeta, dat$gamma)
#gap.plot(x, y, gap=c(200, 3200), type='l',
curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, 3300), #yMax),
        ylab="Yield",
        xlab="B",
        main="Estimated Yield Curves For Poorly Specified BH", #TeX(sprintf("$F_{MSY}/M$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk), round(zeta, binTrk))),
        cex.lab = 1.5,
        cex.main= 1.4
)
polygon(c(grd, rev(grd)), chop(c(qYield[,1], rev(qYield[,3]))),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M), lwd=3, add=T, col=cols[1])
#rug(PBar(dat$alpha, dat$beta, dat$gamma, FMsy(dat$alpha, dat$gamma, dat$M), dat$M), lwd=3)
#rug(PBar(fit$alpha, fit$beta, fit$gamma, FMsy(fit$alpha, fit$gamma, fit$M), fit$M), lwd=3, col='red')
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#legend("topright", legend=c("Schnute Truth", "Low Contrast BH Fit", "High Contrast BH Fit"), col=c("black", cols[1:2]), lwd=3)
dev.off()

#RELATIVE

#
lns = mapply(function(a, b){yield(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qYield = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
lnsS = mapply(function(a, b){yield(grd, a, b, fitS$gamma)}, samS[,1], samS[,2])
qYieldS = rowQuantiles(lnsS, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)

#
png(sprintf("yeildRelCurveCompare%sPrettyX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
x = seq(0, xMax, length.out=1000)
y = yield(x, dat$lalpha, dat$lbeta, dat$gamma)
#gap.plot(x, y, gap=c(200, 3200), type='l',
curve(yield(x*P0, dat$lalpha, dat$lbeta, dat$gamma), 0, 1, n=1000,
        lwd=3,
        ylim=c(0, 3300), #yMax),
        ylab="Yield",
        xlab=TeX("$B/B_0$"),
        main="Estimated Yield Curves For Poorly Specified BH", #TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
P0Fit = PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M)
polygon(c(grd/P0Fit, rev(grd/P0Fit)), chop(c(qYield[,1], rev(qYield[,3]))),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
P0FitS = PBar(fitS$alpha, fitS$beta, fitS$gamma, 0, fitS$M)
polygon(c(grd/P0FitS, rev(grd/P0FitS)), chop(c(qYieldS[,1], rev(qYieldS[,3]))),
        col=adjustcolor(cols[2], alpha.f=0.2),
        border=F
)
curve(yield(x*PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M), fit$lalpha, fit$lbeta, fit$gamma), 0, 1, lwd=3, add=T, col=cols[1])
curve(yield(x*PBar(fitS$alpha, fitS$beta, fitS$gamma, 0, fitS$M), fitS$lalpha, fitS$lbeta, fitS$gamma), 0, 1, lwd=3, add=T, col=cols[2])
#rug(PBar(dat$alpha, dat$beta, dat$gamma, FMsy(dat$alpha, dat$gamma, dat$M), dat$M), lwd=3)
#rug(PBar(fit$alpha, fit$beta, fit$gamma, FMsy(fit$alpha, fit$gamma, fit$M), fit$M), lwd=3, col='red')
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
legend("topright", legend=c("Schnute Truth", "Low Contrast BH Fit", "High Contrast BH Fit"), col=c("black", cols[1:2]), lwd=3)
dev.off()

#
png(sprintf("rR0Compare%sPrettyX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma)/SRR(P0, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000, #yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, 2.5), #3300), #yMax),
        ylab="Production",
        xlab="B",
        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma)/SRR(PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M), fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(PBar(dat$alpha, dat$beta, dat$gamma, FMsy(dat$alpha, dat$gamma, dat$M), dat$M), lwd=3)
rug(PBar(fit$alpha, fit$beta, fit$gamma, FMsy(fit$alpha, fit$gamma, fit$M), fit$M), lwd=3, col='red')
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#######################################################################################################################
#
#MID, MID HI
#
#######################################################################################################################


##
#xi = 3.5
#zeta = 0.2
xiMid = 2
zetaMidHi = 0.45
 
#minimize norm
norms = sqrt((xiMid-xis)^2 + (zetaMidHi-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]
##
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zeta, binTrk))
#
fileDat = names(who)                            #sprintf('%s/datGen_xi%s_zeta%s.rda', dir, round(xi, binTrk), round(zet
fileFit = gsub("datGen", "fit", fileDat)

#
dat = readRDS(fileDat)
fit = readRDS(fileFit)
fitS = readRDS('fitSingle_xi1.93_zeta0.45.rda') #fitSingle_xi3.48_zeta0.48.rda')

#
xMax = max(dat$N0, fit$N0, fitS$N0)
grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
CS = fitS$rsCov[who, who]
mS = c(fitS$lalpha, fitS$lbeta)
samS = rmvnorm(MM, mS, CS)

#
#grd = 0:xMax
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)

#
png(sprintf("yeildCurveCompare%sPrettyX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
x = seq(0, xMax, length.out=1000)
y = yield(x, dat$lalpha, dat$lbeta, dat$gamma)
#gap.plot(x, y, gap=c(200, 3200), type='l',
curve(yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, 3300), #yMax),
        ylab="Yield",
        xlab="B",
        main="Estimated Yield Curves For Poorly Specified BH", #TeX(sprintf("$F_{MSY}/M$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk), round(zeta, binTrk))),
        cex.lab = 1.5,
        cex.main= 1.4
)
polygon(c(grd, rev(grd)), chop(c(qYield[,1], rev(qYield[,3]))),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(yield(x, fit$lalpha, fit$lbeta, fit$gamma), 0, PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M), lwd=3, add=T, col=cols[1])
#rug(PBar(dat$alpha, dat$beta, dat$gamma, FMsy(dat$alpha, dat$gamma, dat$M), dat$M), lwd=3)
#rug(PBar(fit$alpha, fit$beta, fit$gamma, FMsy(fit$alpha, fit$gamma, fit$M), fit$M), lwd=3, col='red')
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
#legend("topright", legend=c("Schnute Truth", "Low Contrast BH Fit", "High Contrast BH Fit"), col=c("black", cols[1:2]), lwd=3)
dev.off()

#RELATIVE

#
lns = mapply(function(a, b){yield(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qYield = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
lnsS = mapply(function(a, b){yield(grd, a, b, fitS$gamma)}, samS[,1], samS[,2])
qYieldS = rowQuantiles(lnsS, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))
##
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qSRR, na.rm=T)

#NOTE: this one
png(sprintf("yeildRelCurveCompare%sPrettyX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
x = seq(0, xMax, length.out=1000)
y = yield(x, dat$lalpha, dat$lbeta, dat$gamma)
#gap.plot(x, y, gap=c(200, 3200), type='l',
curve(yield(x*P0, dat$lalpha, dat$lbeta, dat$gamma), 0, 1, n=1000,
        lwd=3,
        ylim=c(0, 3300), #yMax),
        ylab="Yield",
        xlab=TeX("$B/B_0$"),
        main="Estimated Yield Curves For Poorly Specified BH", #TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
P0Fit = PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M)
polygon(c(grd/P0Fit, rev(grd/P0Fit)), chop(c(qYield[,1], rev(qYield[,3]))),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
P0FitS = PBar(fitS$alpha, fitS$beta, fitS$gamma, 0, fitS$M)
polygon(c(grd/P0FitS, rev(grd/P0FitS)), chop(c(qYieldS[,1], rev(qYieldS[,3]))),
        col=adjustcolor(cols[2], alpha.f=0.2),
        border=F
)
curve(yield(x*PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M), fit$lalpha, fit$lbeta, fit$gamma), 0, 1, lwd=3, add=T, col=cols[1])
curve(yield(x*PBar(fitS$alpha, fitS$beta, fitS$gamma, 0, fitS$M), fitS$lalpha, fitS$lbeta, fitS$gamma), 0, 1, lwd=3, add=T, col=cols[2])
#rug(PBar(dat$alpha, dat$beta, dat$gamma, FMsy(dat$alpha, dat$gamma, dat$M), dat$M), lwd=3)
#rug(PBar(fit$alpha, fit$beta, fit$gamma, FMsy(fit$alpha, fit$gamma, fit$M), fit$M), lwd=3, col='red')
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
legend("topright", legend=c("Schnute Truth", "Low Contrast BH Fit", "High Contrast BH Fit"), col=c("black", cols[1:2]), lwd=3)
dev.off()

#
png(sprintf("rR0Compare%sPrettyX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma)/SRR(P0, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000, #yield(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, 2.5), #3300), #yMax),
        ylab="Production",
        xlab="B",
        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
#polygon(c(grd, rev(grd)), c(qYield[,1], rev(qYield[,3])),
#        col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma)/SRR(PBar(fit$alpha, fit$beta, fit$gamma, 0, fit$M), fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(PBar(dat$alpha, dat$beta, dat$gamma, FMsy(dat$alpha, dat$gamma, dat$M), dat$M), lwd=3)
rug(PBar(fit$alpha, fit$beta, fit$gamma, FMsy(fit$alpha, fit$gamma, fit$M), fit$M), lwd=3, col='red')
#rug(dat$N, lwd=3, ticksize=0.05)
##rug(dat$N, line=0.75)
##rug(fit$N, col=cols[1])
legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
dev.off()












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
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)
#
##
#png(sprintf("curveCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#      lwd=3, 
#      ylim=c(0, yMax),
#      ylab="Production",
#      xlab="B",
#      main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#      cex.lab = 1.5,
#    cex.main= 1.5
#)
#polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
#      col=adjustcolor(cols[1], alpha.f=0.2),
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
#png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
#diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#      lwd=3, 
#      ylim=c(0, yMax),
#      ylab="Production",
#      xlab="B",
#      main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk))),
#      cex.lab = 1.5,
#    cex.main= 1.5
#)
#step=-10
#quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
#      step, step*diffDat,
#      scale=100,
#      col=1:length(dat$N) 
#)
#diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
#      step, step*diffFit,
#      scale=100,
#      col=1:length(dat$N) 
#)
#rug(dat$N, line=0.75)
#rug(fit$N, col=cols[1])
#dev.off()

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
#                       colQuantiles(NSam, probs=0.025),
#                       rev(colQuantiles(NSam, probs=0.975))
#               ),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
#
#
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
#                       colQuantiles(depSam, probs=0.025),
#                       rev(colQuantiles(depSam, probs=0.975))
#               ),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##dPost = t(mapply( function(a, b){
##      #
##      tmp$lalpha = a
##      tmp$lbeta  = b
##      tmp$iterate('lsode')
##      #
##      return( tmp$N/tmp$N[1] )
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
#                       colQuantiles(dPost, probs=0.025),
#                       rev(colQuantiles(dPost, probs=0.975))
#               ),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
###
##bPost = t(mapply( function(a, b){
##      #
##      tmp$lalpha = a
##      tmp$lbeta  = b
##      tmp$iterate('lsode')
##      #
##      return( tmp$N )
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
#                       colQuantiles(bPost, probs=0.025),
#                       rev(colQuantiles(bPost, probs=0.975))
#               ),
#                col = makeTransparent('red', alpha=100),
#                border = NA
#)
#legend("topright", legend=c("Schnute Truth", "BH Fit"), col=c("black", 'red'), lwd=3)
#dev.off()
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$al
##writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.3f\t", dat$xi*M, dat$zeta, (fit$lalpha-dat$lalpha)/sqrt(C['lalpha', 'lalp
#writeLines(sprintf("(%1.1f, %1.1f)\t& (%1.2f, %1.2f, %1.1f)\t& (%1.0f, %1.0f, %1.1f)\t& %1.3f\t", dat$xi*M, dat$zeta, dat$alpha, fit$
#br = ff*fit$catch*fit$N
#dbr = dff*dat$catch*dat$N
#lbr = TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk)*M, round(zeta, binTrk)))






