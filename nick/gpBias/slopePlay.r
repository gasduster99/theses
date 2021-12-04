rm(list=ls())

#
library(pracma)
library(mvtnorm)
library(latex2exp)
library(matrixStats)
library(RColorBrewer)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
SRR = function(P, lalpha, lbeta, gamma){
        #
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
}

#
PBar = function(ff, alpha, beta, gamma, M){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }

#
FMsy = function(alpha, gamma, M){
        ##
        #FUpper = P0 #alpha-M
        ##root = uniroot.all(function(ff){  ((gamma-1)/alpha)^(gamma-1) - (M+ff)^(
        #root = uniroot.all(function(ff){ 1 - ((M+ff)*((gamma-1)/alpha))^(1/(gamma
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        #
        root = alpha/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
        return(root)
}

#
diff = function(f, x, h=0.00001){
        (f(x)-f(x-h))/h
}

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
        Fmsy = exp(lalpha)/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
        out = R - Fmsy*C*P
        #
        return( list(out) )
}

#
#MAIN
#

#
mod = "ExpNoQ"
dir = sprintf("./modsPella%s", mod)
M = 0.2
P0 = 10000

#
#xi = 0.5                #2.3    #
#zeta = 0.2      #0.25   #
Fmsy = 0.1
xi = Fmsy/M
zeta = 0.6
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)

#
diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)

#
diffSam = t(mapply(function(la, lb){ diff(function(x){SRR(x, la, lb, fit$gamma)}, fit$N) }, sam[,1], sam[,2]))

##
#hist(diffSam, 
#	xlim=c(min(diffDat[-1], diffSam), max(diffDat[-1], diffSam)),
#	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
#	xlab="SRR Derivative"
#)
#rug(diffDat)
#rug(diffFit, col='red', line=0.75)
#print(diffDat[diffDat==max(diffDat)])
#print(diffDat[diffDat==min(diffDat)])

#
plot(diffDat[-1], diffFit[-1])
lines(c(-1000, 1000), c(-1000, 1000))







