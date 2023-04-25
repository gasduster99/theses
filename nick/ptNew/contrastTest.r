rm(list=ls())

#
library(tgp)
library(lamW)
library(geoR)
library(VGAM)
library(pracma)
library(foreach)
library(parallel)
library(rootSolve)
library(doParallel)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
dPdt = function(t, P, lalpha, lbeta, gamma, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        Fmsy = exp(lalpha)/gamma 					#exp(lalpha)/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta))^(gamma-1))	#exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
        out = R - Fmsy*C*P
        #
        return( list(out) )
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
fContrast = function(con){
	TT = 31 #length(hake)
	tt = 1
	time = tt:TT
	#
	mid = round(TT/2)
	cMax = 2^con	#2 	#4*M
	cMin = 0.2^con	#0.2 	#M/4
	b = log(cMin/cMax)/(tt-mid)
	a = exp(log(cMax)-b*mid)
	rSlope = (cMax-1)/(mid-1)
	rb = 2*rSlope*mid + cMax-rSlope*mid
	FtFmsy = (a*exp(b*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
	#
	FtFmsy = c(FtFmsy, rep(1, 15))	
}

#
#INIT
#

###
###hake  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43,
###catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 1
##TT = 31 #length(hake)
##tt = 1
##time = tt:TT
###
##FtFmsy = rep(1, TT)
##TT = length(FtFmsy)
##time = tt:TT
#TT = 31 #length(hake)
#tt = 1
#time = tt:TT
##
#mid = round(TT/2)
#con = 2
#cMax = 2^con	#2 	#4*M
#cMin = 0.2^con	#0.2 	#M/4
#b = log(cMin/cMax)/(tt-mid)
#a = exp(log(cMax)-b*mid)
#rSlope = (cMax-1)/(mid-1)
#rb = 2*rSlope*mid + cMax-rSlope*mid
#FtFmsy = (a*exp(b*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
##
#FtFmsy = c(FtFmsy, rep(1, 15))

#con in (0, 2)
con = 1
FtFmsyB = fContrast(con)
TT = length(FtFmsyB)
time = 1:TT #tt:TT

##
#plot(FtFmsy)

#
M = 0.2
#
P0 = 10000 #3000

#
#SIM
#

###grid for simulation
##zetaSims = seq(0.4, 0.5, 0.025) #(seq(0.1, 0.9, 0.01)) #rev(seq(0.15, 0.7, 0.01)) #rev(seq(0.
##xiSims =   (seq(3, 4, 0.25)) #(seq(0.5, 4.5, 0.05))  #rev(seq(0.5, 3.5, 0.05))               
#zetaSims = (seq(0.15, 0.7, 0.05))    #0.65 #seq(0.1, 0.9, 0.05) #(seq(0.1, 0.9, 0.01)) #rev(se
#xiSims = rev(seq(0.5, 3.5, 0.05))  #rev(seq(0, 4, 0.5)) #(seq(0.5, 4.5, 0.05))  #rev(seq(0.5, 
##zetaSims = (seq(0.45, 0.65, 0.1))      
##xiSims = rev(seq(0.5, 4.5, 0.05))
zetaMax = 0.7
zetaMin = 0.15
xiMax = 0.8  #3.75
xiMin = 0.05 #0.25
#
n = 100
minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
binTrk = ceiling(abs(log10(minDiff)))

#modsPTFlatT30/fit_xi0.699_zeta0.201.rda
xiBase = 0.699
zetaBase = 0.201
baseFitName = sprintf('./modsPTFlatT30/fit_xi%s_zeta%s.rda', round(xiBase, binTrk), round(zetaBase, binTrk))
#
baseDatName = gsub("fit", "datGen", baseFitName)
baseDat = readRDS(baseDatName)
baseFit = readRDS(baseFitName)

#a place to store data
place = sprintf('./contrastTestsX%sZ%s/', round(xiBase, binTrk), round(zetaBase, binTrk))
odeMethod = "lsode"


#make new design
if( F ){
	#
	cons = seq(0, 2, length.out=46)
	#
	for(c in cons){
		#fishing
		FtFmsy = fContrast(c)
		TT = length(FtFmsy)
		time = 1:TT #tt:TT
		#read in data from design locations
	        datGen = prodModel$new(
	                dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #
	                time=1:TT, catch=FtFmsy,  					#constants
	                alpha=baseDat$alpha, beta=baseDat$beta, gamma=baseDat$gamma,	#parameters
	                lalpha=baseDat$lalpha, lbeta=baseDat$lbeta,      	 	#reparameterize
	                lq=baseDat$lq, lsdo=baseDat$lsdo,          			#nuisance parameters
	                xi=baseDat$xi, zeta=baseDat$zeta, con=c		     	#other incidentals to carry along
	        )
	        datGen$iterate(odeMethod)
		datName = sprintf('%s/datGen_xi%s_zeta%s_con%s.rda', place, round(datGen$xi, binTrk), round(datGen$zeta, binTrk), round(c, binTrk))
	        datGen$save( datName )
	}
}else{

#
datFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("datGen*.rda")))

#
registerDoParallel(46)
opts = list(preschedule=F)
foreach(i=(1:length(datFiles)), .options.multicore = opts) %dopar% {
#for(i in 1:length(datFiles)){
	#
        #DATA
        #

	#read in data from design locations
        datGen = readRDS(datFiles[i])
        #datGen$time = 1:TT
        #datGen$catch = FtFmsy
        #datGen$iterate(odeMethod)
        #datGen$save( datFiles[i] )
		
	#
        #SKIP
        #

        #
        fileFit = gsub("datGen", "fit", datFiles[i])
        if( file.exists(fileFit) ){
                #writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
                next
        }

        #
        #DATA
        #

        #make data to fit
        cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
	#datGen$plotMean()
        #points(cpue)
	
	#
        #FIT
        #

        #NOTE: N0Funk does not function correctly with $clone()
        fit = prodModel$new(
                dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #
                time=1:TT, catch=datGen$catch,    #schaefer   #constants
                alpha=datGen$alpha, beta=datGen$beta, gamma=2, 	#parameters
                lalpha=datGen$lalpha, lbeta=datGen$lbeta,  	#reparameterize
                lq=datGen$lq, lsdo=datGen$lsdo, 	#nuisance parameters
                xi=datGen$xi, zeta=datGen$zeta		#other incidentals to carry along
        )
        #fit = readRDS('modsPellaFineQFixRFixP010000/fit_xi4_zeta0.35.rda')
        fit$iterate(odeMethod)

	#optimization
        optAns = fit$optimize(cpue,
                c('lsdo', 'lalpha'), #'lq'),
                lower   = c(log(0.001), log(M)), #log(1e-7)),
                upper   = c(log(1), log(100)), #log(1e-2)),
                gaBoost = list(run=10, parallel=FALSE, popSize=10^3),
                persistFor = 5,
                fitQ    = F
        )
        optAns = fit$optimize(cpue,
                c('lsdo', 'lalpha', 'lbeta'),
                lower   = c(log(0.001), log(0), log(P0/10)),    #log(10^3)), 
                upper   = c(log(1), log(100), log(1.5*P0)),     #log(10^4)),
                gaBoost = list(run=100, parallel=FALSE, popSize=10^3),
                persistFor = 5,
                fitQ    = F
        )
	#get hessian if possible
        tryCatch({
                optAns = fit$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta'),
                        lower   = c(log(0.001), log(0), log(P0/10)),
                        upper   = c(log(1), log(100), log(1.5*P0)),
                        cov     = T,
                        fitQ    = F
                        #c('lsdo', 'alpha', 'beta', 'lq'),                                                      
                        #lower   = c(log(0.001), eps(), eps(), log(1e-7)),
                        #upper   = c(log(1), 10, 2, log(1e-2)),
                )
        }, error=function(err){
                writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
                optAns = fit$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta'),
                        lower   = c(log(0.001), log(0), log(P0/10)),
                        upper   = c(log(1), log(100), log(1.5*P0)),
                        cov     = F,
                        fitQ    = F
                        #c('lsdo', 'alpha', 'beta', 'lq'),                                                      
                        #lower   = c(log(0.001), eps(), eps(), log(1e-7)),
                        #upper   = c(log(1), 10, 2, log(1e-2)), 
                )
        })
	#convenience
        fit$alpha = exp(fit$lalpha)
        fit$beta  = exp(fit$lbeta)
        fit$sdo   = exp(fit$lsdo)
        fit$q     = exp(fit$lq)
        #save
        fit$save( fileFit )
}


}


