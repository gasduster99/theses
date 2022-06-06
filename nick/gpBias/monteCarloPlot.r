rm(list=ls())

#
library(mvtnorm)
library(rootSolve)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

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
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) }

#
getAlpha = function(gamma, ff, M){
        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}
a = Vectorize(getAlpha, "gamma")

#
getBeta = function(alpha, gamma, M, B0){
        (1-(M/alpha)^gamma)/B0/gamma
}
b = Vectorize(getBeta, c("alpha", "gamma"))

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
getFits = function(dir, xi, zeta, binTrk=2){
	#
	l = unlist(strsplit(gsub('.rda', '', names(xi)), '_'))
	whoX = grep('xi', l)
	whoZ = grep('zeta', l)
	xiBin = gsub('xi', '', l[whoX])
	zetaBin = gsub('zeta', '', l[whoZ])
	#
	globber = sprintf("fit_xi%s_zeta%s_*.rda", xiBin, zetaBin)
        #
	sprintf( "%s%s", dir, list.files(path=dir, pattern=glob2rx(globber)) )
}

#
getData = function(dir, xi, zeta){
        #
	fitFiles = getFits(dir, xi, zeta)
       	#
        i = 1
        D = data.frame(xiBin=double(), zetaBin=double(), xiHat=double(), zetaHat=double(), lF=double(), lFV=double(), lK=double(), lKV=double(), likelihood=double(), stringsAsFactors=F)
        for(f in fitFiles){
                #
                fit = readRDS(f) 
		fitLike = prodModel$new(
		        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},        #model
		        time=fit$time, catch=fit$catch, M=fit$M,      	    #constants
		        alpha=fit$alpha, beta=fit$beta, gamma=-1,           #parameters
		        lalpha=fit$lalpha, lbeta=fit$lbeta,                 #reparameterize
		        lq=fit$lq, lsdo=fit$lsdo,                           #nuisance parameters
		        xi=fit$xi, zeta=fit$zeta                            #other incidentals to carry along
		)
		fitLike$iterate('lsode')
		fitLike$iterate('lsode')

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

                ##
                #md = optimize(myDist, c(0, max(xiRange)), xi=dat$xi, zeta=dat$zeta)$objective

                #
                if(length(fit$rsCov)==0){ lfv=0; lbv=0 }else{
                        #
                        lV = getlV(fit)
                        #
                        lfv = lV$lFV
                        lbv = lV$lKV
                }

                #
                D[i,] = c(xiBin, zetaBin, xiHat, zetaHat, log(Fs), lfv, log(BZero), lbv, fitLike$like(fit$cpue))
                i = i+1
        }
        #
        return(D)
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
#
#

#
M = 0.2
P0 = 10000
#
mod = "HHardFlatT30N150Wide"
dirOut = sprintf("./monteCarlo%s/", mod)
dirIn = sprintf("./modsSchnute%s/", mod)

#
datFiles = sprintf("%s%s", dirIn, list.files(path=dirIn, pattern=glob2rx("datGen*.rda")))
#
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

##min, min
#xiPoint = 0.5 		#(3.5-0.5)/2 + 0.5
#zetaPoint = 0.2 	#(0.7-0.2)/2 + 0.2
##
#norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
#who   = which(min(norms)==norms)
#xiMinMin = xis[who]
#zetaMinMin  = zetas[who]

#max, mid
xiPoint = 3.46
zetaPoint = 0.4
#
norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
who   = which(min(norms)==norms)
xiMaxMid = xis[who]
zetaMaxMid  = zetas[who]

##min, mid
#xiPoint = 0.5
#zetaPoint = 0.593
##
#norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
#who   = which(min(norms)==norms)
#xiMinMid = xis[who]
#zetaMinMid  = zetas[who]
#
##max, max
#xiPoint = 3.5
#zetaPoint = 0.7
##
#norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
#who   = which(min(norms)==norms)
#xiMaxMax = xis[who]
#zetaMaxMax  = zetas[who]
#
##mid, mid
#xiPoint = (3.5-0.5)/2 + 0.5
#zetaPoint = (0.7-0.2)/2 + 0.2
##
#norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
#who   = which(min(norms)==norms)
#xiMidMid = xis[who]
#zetaMidMid  = zetas[who]
#
##
#png(sprintf('%sDesignMod.png', mod))
#plot(xis, zetas, pch=20)
##abline(v=xiL, lty=3)
##abline(h=zetaL, lty=3)
#curve(1/(x+2), 0, 4, add=T)
#points(xiMidMid, zetaMidMid, pch=4, col=6, cex=1.5)
#points(xiMaxMax, zetaMaxMax, pch=4, col=2, cex=1.5)
#points(xiMinMid, zetaMinMid, pch=4, col=3, cex=1.5)
#points(xiMaxMid, zetaMaxMid, pch=4, col=4, cex=1.5)
#points(xiMinMin, zetaMinMin, pch=4, col=5, cex=1.5)
#dev.off()


##MidMid
#fitNames = getFits(dirOut, xiMidMid, zetaMidMid)
#fits = sapply(fitNames, function(fn){
#	readRDS(fn)
#})
#dat  = getData(dirOut, xiMidMid, zetaMidMid)
##
#png(sprintf('%sXiHatMidMid.png', mod))
#hist(dat$xiHat, col=6, xlim=c(min(dat$xiHat, xiMidMid), max(dat$xiHat, xiMidMid)))
#abline(v=xiMidMid)
#dev.off()
##
#png(sprintf('%sZetaHatMidMid.png', mod))
#hist(dat$zetaHat, col=6, xlim=c(min(dat$zetaHat, zetaMidMid), max(dat$zetaHat, zetaMidMid)))
#abline(v=zetaMidMid)
#dev.off()
##
#png(sprintf('%sFmsyHatMidMid.png', mod))
#hist(exp(dat$lF), col=6, xlim=c(min(exp(dat$lF), xiMidMid*M), max(dat$lF, xiMidMid*M)))
#abline(v=xiMidMid*M)
#dev.off()
##
#png(sprintf('%sKHatMidMid.png', mod))
#hist(exp(dat$lK), col=6, xlim=c(min(exp(dat$lK), P0), max(exp(dat$lK), P0)))
#abline(v=P0)
#dev.off()

##
#i = 1
#fit = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},	#model
#        time=fits[[i]]$time, catch=fits[[i]]$catch, M=fits[[i]]$M,  	#constants
#        alpha=fits[[i]]$alpha, beta=fits[[i]]$beta, gamma=-1, 		#parameters
#        lalpha=fits[[i]]$lalpha, lbeta=fits[[i]]$lbeta,       		#reparameterize
#        lq=fits[[i]]$lq, lsdo=fits[[i]]$lsdo,    			#nuisance parameters
#        xi=fits[[i]]$xi, zeta=fits[[i]]$zeta                  		#other incidentals to carry along
#)
#fit$iterate('lsode')

##MaxMax
#fitNames = getFits(dirOut, xiMaxMax, zetaMaxMax)
#fits = sapply(fitNames, function(fn){
#	readRDS(fn)
#})
#dat  = getData(dirOut, xiMaxMax, zetaMaxMax)
##
#png(sprintf('%sXiHatMaxMax.png', mod))
#hist(dat$xiHat, col=2, xlim=c(min(dat$xiHat, xiMaxMax), max(dat$xiHat, xiMaxMax)))
#abline(v=xiMaxMax)
#dev.off()
##
#png(sprintf('%sZetaHatMaxMax.png', mod))
#hist(dat$zetaHat, col=2, xlim=c(min(dat$zetaHat, zetaMaxMax), max(dat$zetaHat, zetaMaxMax)))
#abline(v=zetaMaxMax)
#dev.off()
##
#png(sprintf('%sFmsyHatMaxMax.png', mod))
#hist(exp(dat$lF), col=2, xlim=c(min(exp(dat$lF), xiMaxMax*M), max(exp(dat$lF), xiMaxMax*M)))
#abline(v=xiMaxMax*M)
#dev.off()
##
#png(sprintf('%sKHatMaxMax.png', mod))
#hist(exp(dat$lK), col=2, xlim=c(min(exp(dat$lK), P0), max(exp(dat$lK), P0)))
#abline(v=P0)
#dev.off()
#
##
#fitNames = getFits(dirOut, xiMinMid, zetaMinMid)
#fits = sapply(fitNames, function(fn){
#	readRDS(fn)
#})
#dat  = getData(dirOut, xiMinMid, zetaMinMid)
##
#png(sprintf('%sXiHatMinMid.png', mod))
#hist(dat$xiHat, col=3, xlim=c(min(dat$xiHat, xiMinMid), max(dat$xiHat, xiMinMid)))
#abline(v=xiMinMid)
#dev.off()
##
#png(sprintf('%sZetaHatMinMid.png', mod))
#hist(dat$zetaHat, col=3, xlim=c(min(dat$zetaHat, zetaMinMid), max(dat$zetaHat, zetaMinMid)))
#abline(v=zetaMinMid)
#dev.off()
##
#png(sprintf('%sFmsyHatMinMid.png', mod))
#hist(exp(dat$lF), col=3, xlim=c(min(exp(dat$lF), xiMinMid*M), max(exp(dat$lF), xiMinMid*M)))
#abline(v=xiMinMid*M)
#dev.off()
##
#png(sprintf('%sKHatMinMid.png', mod))
#hist(exp(dat$lK), col=3, xlim=c(min(exp(dat$lK), P0), max(exp(dat$lK), P0)))
#abline(v=P0)
#dev.off()

#
fitNames = getFits(dirOut, xiMaxMid, zetaMaxMid)
fits = sapply(fitNames, function(fn){
        readRDS(fn)
})
dat  = getData(dirOut, xiMaxMid, zetaMaxMid)
#
png(sprintf('%sXiHatMaxMid.png', mod))
hist(dat$xiHat, col=4, xlim=c(min(dat$xiHat, xiMaxMid), max(dat$xiHat, xiMaxMid)))
abline(v=xiMaxMid)
dev.off()
#
png(sprintf('%sZetaHatMaxMid.png', mod))
hist(dat$zetaHat, col=4, xlim=c(min(dat$zetaHat, zetaMaxMid), max(dat$zetaHat, zetaMaxMid)))
abline(v=zetaMaxMid)
dev.off()
#
png(sprintf('%sFmsyHatMaxMid.png', mod))
hist(exp(dat$lF), col=4, xlim=c(min(exp(dat$lF), xiMaxMid*M), max(exp(dat$lF), xiMaxMid*M)))
abline(v=xiMaxMid*M)
dev.off()
#
png(sprintf('%sKHatMaxMid.png', mod))
hist(exp(dat$lK), col=4, xlim=c(min(exp(dat$lK), P0), max(exp(dat$lK), P0)))
abline(v=P0)
dev.off()
#
png(sprintf('%sLikeHist.png', mod))
hist(dat$likelihood, main="Likelihood at Max", xlab='', col=4)
dev.off()

##
#fitNames = getFits(dirOut, xiMinMin, zetaMinMin)
#fits = sapply(fitNames, function(fn){
#        readRDS(fn)
#})
#dat  = getData(dirOut, xiMinMin, zetaMinMin)
##
#png(sprintf('%sXiHatMinMin.png', mod))
#hist(dat$xiHat, col=5, xlim=c(min(dat$xiHat, xiMinMin), max(dat$xiHat, xiMinMin)))
#abline(v=xiMinMin)
#dev.off()
##
#png(sprintf('%sZetaHatMinMin.png', mod))
#hist(dat$zetaHat, col=5, xlim=c(min(dat$zetaHat, zetaMinMin), max(dat$zetaHat, zetaMinMin)))
#abline(v=zetaMinMin)
#dev.off()
##
#png(sprintf('%sFmsyHatMinMin.png', mod))
#hist(exp(dat$lF), col=5, xlim=c(min(exp(dat$lF), xiMinMin*M), max(exp(dat$lF), xiMinMin*M)))
#abline(v=xiMinMin*M)
#dev.off()
##
#png(sprintf('%sKHatMinMin.png', mod))
#hist(exp(dat$lK), col=5, xlim=c(min(exp(dat$lK), P0), max(exp(dat$lK), P0)))
#abline(v=P0)
#dev.off()





