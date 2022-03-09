rm(list=ls())

#
#library(geoR)
library(sp)
library(gstat)
#library(RandomFields)
#
library(VGAM)
library(boot)
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
getPar = function(xi, zeta, M){
        #
        gamma = 1/zeta
        odd = (1-zeta)/zeta
        alpha = xi*M*odd*(1-zeta)^(-odd)
        #
        return(c(alpha, gamma))
}

#
myDist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/2)^2) }

#
getData = function(dir, xiSims, zetaSims){
        #dir    : a directory containing data
        #xiSims : xis of simulated data
        #zetaSims: zetas of simulated data
        #
        #Seed values are the initiating values to launch the inversion
        #Inv values are the numerical inversions actually found for the 3 parameter curve at each seed value
        #BH values are the MLE fits under the BH model.

        #
        di = 1
        D = data.frame(xiSeed=double(), zetaSeed=double(), xiInv=double(), zetaInv=double(), xiHat=double(), zetaHat=double(), minDist=double(), lF=double(), lFV=double(), lK=double(), lKV=double(), stringsAsFactors=F)
        for(i in 1:length(zetaSims)){
                for(j in 1:length(xiSims)){
                        #
                        fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
                        fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])

                        #
                        if( file.exists(fileDat) & file.exists(fileFit) ){
                                #
                                dat = readRDS(fileDat)
                                fit = readRDS(fileFit)

                                # length(dat$time)!=length(dat$N) |
                                if( length(fit$time)!=length(fit$N) ){ next }

                                #the inversion actually found
                                FInv = FMsy(dat$alpha, dat$gamma, M)
                                xiInv = FInv/M
                                zetaInv = PBar(FInv, dat$alpha, dat$beta, dat$gamma, M)/PBar(0, dat$alpha, dat$beta, dat$gamma, M)
                                #the bh fit found
                                Fs = FMsy(fit$alpha, 2, M)
                                xiHat   = Fs/M
                                zetaHat = PBar(Fs, fit$alpha, fit$beta, fit$gamma, M)/PBar(0, fit$alpha, fit$beta, fit$gamma, M) #xiHat/(2*xiHat+1) #
                                md = optimize(myDist, c(0, max(xiSims)), xi=dat$xi, zeta=dat$zeta)$objective

                                #NOTE: replace with a propper observation uncertainty
                                if(length(fit$rsCov)==0){ v=0; lbv=0 }else{ v=getlFV(fit)$lFV; lbv=fit$rsCov['lbeta', 'lbeta'] }
                                #print( c(dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), v) )
                                D[di,] = c(dat$xi, dat$zeta, xiInv, zetaInv, xiHat, zetaHat, md, log(Fs), v, fit$lbeta, lbv)
                                #print(dim(D))
                                di = di+1
                        }
                }
        }
        #
        return(D)
}

#
getlFV = function(fit, MM=10^4, samples=F){
        #
        who = c('lalpha')
        C = fit$rsCov[who, who]
        m = c(fit$lalpha)
        sam = rnorm(MM, m, sqrt(C)) #rmvnorm(M, m, C)
        #
        als = exp(sam)
        lFs = sapply(als, function(a){log(FMsy(a, 2, M))})
        #
        out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T))
        if( samples ){ out$lFSamples=lFs }
        #
        return(out)
}

#
#TRY
#

#
P0 = 10000
mod = "ExpT45" #"ExpK1T60" #"ExpK" #"Exp1KT60" 
dir = sprintf("./modsPella%s", mod)

#
xiRes = 0.5
zetaTop = 0.6 #0.7
zetaBot = 0.2 #0.1
xiBot = 0.5
xiTop = 3.5
#
zetaSims = seq(0.15, 0.7, 0.1) #0.01 #seq(0.1, 0.8, 0.05)       #rev(seq(
xiSims =   seq(0, 4.5, 0.5) #0.05 #seq(0.5, 3.5, 0.25)          #seq(0.5,
#zetaSims = seq(0.15, 0.9, 0.025)
#xiSims =   seq(0, 4.5, 0.25) #0.

#
M = 0.2
#time: 70
Dall = getData(dir, xiSims, zetaSims)
D = Dall[Dall$lF<4,]

#
lf = as.data.frame(D[,c('xiInv', 'zetaInv', 'lF')])
lk = as.data.frame(D[,c('xiInv', 'zetaInv', 'lK')])
#
coordinates(lf) = ~xiInv+zetaInv 
coordinates(lk) = ~xiInv+zetaInv 

#
g = gstat(NULL, id='lF', form=lF~1, data=lf)
g = gstat(g, id='lK', form=lK~1, data=lk)
vCloud = variogram(g, cutoff=10, cloud=T)
vCross = vCloud[vCloud$id=='lF.lK',]
vCross = vCross[vCross$gamma<1 & vCross$gamma>-1,]

#
png(sprintf('crossCov%s.png', mod))
plot(vCross$dist, vCross$gamma, ylim=c(-1, 1), main="Cross Covariogram", xlab="Distance", ylab="gamma")
low = lowess(vCross$dist, vCross$gamma)
lines(low, col='red')
dev.off()

#line = lm(vCross$gamma~vCross$dist)
#abline(a=line$coef[1], b=line$coef[2], col='red')


#plot(vCross)


##
#dat = as.matrix(D[,c('xiInv', 'zetaInv', 'lF', 'lK')])
#di = as.matrix(dist(D[,c('xiInv', 'zetaInv')]))
#ds = sort(unique(di))
##
#gr = data.frame(c1=D[,'xiInv'], c2=D[,'zetaInv'])
#gridded(gr) = ~c1+c2
##
#


#which(as.matrix(di)==ds[1], arr.ind=T)
##for(d in ds){
##	which(di==d, arr.ind=T)	
##}


####
###dat = as.matrix(D[,c('xiInv', 'zetaInv', 'lF', 'lK')])
###dat = cbind(dat, 1)
##gDF = as.geodata(D[,c('xiInv', 'zetaInv', 'lF')])
###VF  = variog(gDF, op='cloud') 
##VF = variog(gDF, uvec=10)
###plot(VF$u, VF$v)#, ylim=c(0, 1))
##
###
##dev.new()
##gDK = as.geodata(D[,c('xiInv', 'zetaInv', 'lK')])
##VK  = variog(gDK, uvec=10) #max.dist=30
##plot(VK$u, VK$v, ylim=c(0, 1))
#
##
#model <- RMbiwm(A=matrix(c(1,1,1,2), nc=2),
#                nudiag=c(0.5,2), s=c(3, 1, 2), c=c(1, 0, 1))
#x <- seq(0, 20, 1)
#dta <- RFsimulate(model, x, x, n=1)
#d = RFspDataFrame2conventional(dta)
#
###
##dat = 
##
#
##
#
#
###
##crossVario = RFvariogram(data=dat) #, phi=2)
##plot(crossVario)

