rm(list=ls())

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
                D[i,] = c(xiBin, zetaBin, dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), lfv, log(BZero), lbv)
                i = i+1
        }
        #
        return(D)
}

#
getDataSeq = function(dir, xiRange, zetaRange, it=1){
        #

        #
        fitFiles = sprintf( "%s%s", dir, list.files(path=dir, pattern=glob2rx("fit*.rda")) )
	fitFNames = sapply(strsplit(fitFiles, '/'), function(x){x[length(x)]})
	#
	allApps = read.csv(sprintf("%s/appends.csv", dir))
	allAppFNames = sapply(strsplit(allApps$modName, '/'), function(x){x[length(x)]})
	#
	fitFiles = fitFiles[!(fitFNames%in%allAppFNames)]
	fitFiles = c(fitFiles, allApps[allApps$batch<=it,'modName'])	
	#fitFiles = allApps[allApps$batch=it,'modName']
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
                md = optimize(myDist, c(0, max(xiRange)), xi=dat$xi, zeta=dat$zeta)$objective

                #
		#print(fit$rsCov)
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
#DATA
#

#
P0 = 10000
mod = "HHardFlatT30N150WWideN56" #"HHardExpT45N150M0.1Wide" # #"ExpT45N150M0.3Wide" #"HHardFlatT30N150WWideExtra" # #"Exp
place = sprintf("./modsSchnute%s/", mod)

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
sd2lF = c()
obslF = c()
sd2lK = c()
obslK = c()
for(it in 0:56){
	#
	#D = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
	D = getDataSeq(place, c(xiBot, xiTop), c(zetaBot, 0.7), it=it)
	print(dim(D))
	D = D[D$lFV>0 & D$lKV>0,]
	print(dim(D))
	#D = D[c(rep(T, 5), F),]
	#D = D[seq(1, nrow(D), 2),]
	#D = Dall[Dall$lF<4,]
	#plot(D[,1], D[,2], pch=20)
	#points(D[,3], D[,4])
	
	#
	#GP INTERPOLATION
	#
	
	#First fit the lF Model and make lF related Predictions
	
	#pick a polynomial mean function
	lFy = D$lF
	lFV = diag(D$lFV)
	lFX = cbind(1, D$xiSeed, D$zetaSeed)
	
	#
	xAug = seq(0.5, 4, 0.25) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
	aug = cbind(rep(1, length(xAug)), xAug, 1/(xAug+2))
	lFX = rbind(lFX, aug)
	lFy = c(lFy, log(xAug*M))
	lFV = diag(c(D$lFV, rep(mean(D$lFV), length(xAug))))
	
	#
	lFaxes = lFX[,2:3]
	registerData(lFy, lFX, lFaxes, lFV)
	par = c(l1=0.5, l2=0.5, s2=0.5)
	lFFit = gpMAP(par, hessian=F, psiSample=F, lower=c(rep(eps(), 3))) #lower=c(0.3, rep(eps(), 2)))	
	sd2lF = c(sd2lF, lFFit$psi['s2'])
	obslF = c(obslF, mean(D$lFV))
	
	#
	lKy = D$lK
	lKV = diag(c(D$lKV))
	lKX = cbind(1, D$xiSeed, D$zetaSeed)
	
	##
	#xAug = seq(0.75, 4, 0.5) #xAug = c(xAug, seq(7/8, 3, xiRes)) #seq(7/8, 4.5, xiRes)
	#aug = cbind(rep(1, length(xAug)), xAug, 1/(xAug+2))
	#lKX = rbind(lKX, aug)
	#lKy = c(lKy, log(P0))
	#lKV = diag(c(D$lKV, rep(0, length(xAug))))
	
	#
	lKaxes = lKX[,2:3]
	#
	registerData(lKy, lKX, lKaxes, lKV)
	par = c(l1=0.5, l2=0.5, s2=0.5)
	lKFit = gpMAP(par, hessian=F, psiSample=F, lower=c(eps(), rep(eps(), 2)))	
	sd2lK = c(sd2lK, lKFit$psi['s2'])
	obslK = c(obslK, mean(D$lKV))	

	writeLines(sprintf("\nIteration: %s n: %s\n", it, nrow(D)))
	#writeLines("lF Fit:")
	#print(lFFit)
	#print(mean(D$lFV))
	#writeLines('')
	#writeLines("lK Fit:")
	#print(lKFit)
	#print(mean(D$lKV))
}





