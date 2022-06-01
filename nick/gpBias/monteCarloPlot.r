rm(list=ls())

#
#FUNCTIONS
#

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
	globber = sprintf("fit_xi%s_zeta%s_*.rda", round(xi, binTrk), round(zeta, binTrk))
        #
	sprintf( "%s%s", dir, list.files(path=dir, pattern=glob2rx(globber)) )
}

#
getData = function(dir, xi, zeta){
        #
	fitFiles = getFits(dir, xi, zeta)
       	#
        i = 1
        D = data.frame(xiBin=double(), zetaBin=double(), xiHat=double(), zetaHat=double(), lF=double(), lFV=double(), lK=double(), lKV=double(), stringsAsFactors=F)
        for(f in fitFiles){
                #
                fit = readRDS(f) 

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
                D[i,] = c(xiBin, zetaBin, xiHat, zetaHat, log(Fs), lfv, log(BZero), lbv)
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

#
mod = "FlatT30N150Wide"
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

#
xiPoint = (3.5-0.5)/2 + 0.5
zetaPoint = (0.7-0.2)/2 + 0.2
#minimize (max, max) norm
norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]

#
fits = getFits(dirOut, xi, zeta)
dat  = getData(dirOut, xi, zeta)
