rm(list=ls())

#
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

##
#getGamma = function(zeta, ff, M){
#       uniroot.all(function(gamma){z(gamma, ff, M)-zeta}, c(-500, 100))
#       #uniroot.all(function(gamma){ ((ff/(M+ff))*(1-zeta)/(zeta*gamma)+1)^gamma - (M+ff)/M }, c(
#       #strongRoot()
#}

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
#INIT
#

#
#hake  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63, 0.64)
#catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = 31 #length(hake)
tt = 1
time = tt:TT
#
mid = round(TT/2)
cMax = 2 #4*M
cMin = 0.2 #M/4
bb = log(cMin/cMax)/(tt-mid)
aa = exp(log(cMax)-bb*mid)
rSlope = (cMax-1)/(mid-1)
rb = 2*rSlope*mid + cMax-rSlope*mid
FtFmsy = (aa*exp(bb*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
#
FtFmsy = c(FtFmsy, rep(1, 15))
TT = length(FtFmsy)
time = tt:TT
#
M = 0.2
#
P0 = 10000 #3000

#
#SIM
#

#a place to store data
place = './modsSchnuteExpT45N150/'
odeMethod = "lsode"

#
datFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("datGen*.rda")))

#
registerDoParallel(46) 
opts = list(preschedule=F)
foreach(i=rev(1:length(datFiles)), .options.multicore = opts) %dopar% {
#for(i in 1:length(datFiles)){
	#
	#DATA
	#
	
	#read in data from design locations
	datGen = readRDS(datFiles[i])
	datGen$time = 1:TT
	datGen$catch = FtFmsy
	datGen$iterate(odeMethod)
	datGen$save( datFiles[i] )

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
	#writeLines(sprintf("\nPID: %s, Xi: %s, Zeta: %s -> (%f, %f, %f) \n", Sys.getpid(), xiSims[j], zetaSims[i], alpha, beta, gamma))

	#
	#FIT
	#

	#NOTE: N0Funk does not function correctly with $clone()
	fit = prodModel$new(
	        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)}, #model
	        time=1:TT, catch=FtFmsy, M=M,				     #BH	#constants
	        alpha=datGen$alpha, beta=getBeta(datGen$alpha, -1, M, P0), gamma=-1, 	#parameters
		lalpha=datGen$lalpha, lbeta=log(getBeta(datGen$alpha, -1, M, P0)), 	#reparameterize
	        lq=log(0.00049), lsdo=log(0.01160256), 		#nuisance parameters
	        xi=datGen$xi, zeta=datGen$zeta			#other incidentals to carry along
	)
	#fit = readRDS('modsPellaFineQFixRFixP010000/fit_xi4_zeta0.35.rda')
	fit$iterate(odeMethod)
	#optimization	
	optAns = fit$optimize(cpue,
	        c('lsdo', 'lalpha'), 
	        lower   = c(log(0.001), log(M)), 
	        upper   = c(log(1), log(100)),
	        gaBoost = list(run=10, parallel=F, popSize=10^3),
		persistFor = 5,
		fitQ    = F
	)
	#
	optAns = fit$optimize(cpue,
	        c('lsdo', 'lalpha', 'lbeta'), 
	        lower   = c(log(0.001), log(M), -10), 	
	        upper   = c(log(1), log(100), -2), 	#log(getBeta(100, -1, M, P0))
	        gaBoost = list(run=100, parallel=F, popSize=10^3),
		persistFor = 5,
		fitQ    = F
	)
	#get hessian if possible
	tryCatch({
		optAns = fit$optimize(cpue,
			c('lsdo', 'lalpha', 'lbeta'), 
                 	lower   = c(log(0.001), log(M), -10), 
                 	upper   = c(log(1), log(100), -2), 
			cov     = T,
			fitQ    = F
		)
	}, error=function(err){
		writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
		optAns = fit$optimize(cpue,
			c('lsdo', 'lalpha', 'lbeta'), 
			lower   = c(log(0.001), log(M), -10), 
			upper   = c(log(1), log(100), -2), 
			cov     = F,
			fitQ    = F
		)
	})
	#convenience
	fit$alpha = exp(fit$lalpha)
	fit$beta  = exp(fit$lbeta)
	fit$sdo   = exp(fit$lsdo)
	fit$q     = exp(fit$lq) 
	#save
	fit$save( fileFit )
	##plot
	#fit$plotMean(add=T, col='blue')
	#fit$plotBand(col='blue', alpha=50)
}











#
#JUNK
#


	#datGen$printSelf() #catch = 
#	datGen = prodModel$new(
#	        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #model
#	        time=1:TT, catch=FtFmsy, M=M,                    #constants
#	        alpha=alpha, beta=beta, gamma=gamma,            #parameters
#		lalpha=log(alpha), lbeta=log(beta),		#reparameterize
#	        lq=log(0.00049), lsdo=log(0.01160256) #log(0.1160256)		#nuisance parameters
#	        #xi=xi, zeta=zeta                                #other incidentals to carry along
#	)
#	datGen$iterate(odeMethod)
#	
#	#
#        for(j in 1:length(xiSims)){
#		#
#		#SKIP
#		#
#		
#		#
#		fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i])
#		if( file.exists(fileFit) ){ 
#			#writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
#			next 
#		}
#
#		#
#		#INVERT
#		#
#		
#		#
#                Fs = M*xiSims[j]
#                Ps = P0*zetaSims[i]
#
#		##
#		#par = c(alpha, gamma)
#		#extra = data.frame(xi=xiSims[j], zeta=zetaSims[i], M=M, P0=P0)
#		#capture.output( MR <- multiroot(f, par, parms=extra, maxiter=1e6, positive=T), file="/dev/null" )
#		##
#		#par = MR$root
#		#if( !isGood(par, extra) ){
#		#        par = strongRoot(f, par, extra, howGood, lower=c(0, 1), upper=c(5, 5), monitor=F)
#		#}	
#		##
#		#alpha = par[1]
#		#gamma = par[2]
#		#beta  = P0 #getBeta(alpha, gamma, M, P0)
#		##
#		##par values 
#        	#FStar = FMsy(alpha, gamma, M)
#        	#PStar = PBar(FStar, alpha, beta, gamma, M)
#		#PZero = P0 #PBar(0, alpha, beta, gamma, M)
#		
#		#
#		par = getPar(xiSims[j], zetaSims[i], M)
#		#
#		alpha = par[1]
#                gamma = par[2]
#                beta  = P0
#		#
#		FStar = FMsy(par[1], par[2], M) #uniroot.all
#		PStar = PBar(FStar, par[1], beta, par[2], M)
#
#		##
#		#writeLines(sprintf("Xi: %s, Zeta: %s -> (%f, %f, %f) -> (%f, %f) isGood: %d \n", xiSims[j], zetaSims[i], alpha, beta, gamma, FStar/M, PStar/PZero, isGood(par, extra)))
#
#		#
#		#DATA
#		#
#		
#		#update data generator
#		datGen$alpha = alpha;	datGen$lalpha = log(alpha)
#		datGen$beta  = beta;	datGen$lbeta  = log(beta)
#		datGen$gamma = gamma
#		#convience
#		datGen$sdo 	= exp(datGen$lsdo)	
#		datGen$q   	= exp(datGen$lq)
#		datGen$xi	= xiSims[j]
#		datGen$zeta	= zetaSims[i]
#		#iterate
#		datGen$iterate(odeMethod)	
#		
#		##
#		#datGen$printSelf()
#
#		##
#		#plot(c(0,0), ylim=c(0.25, 2), xlim=c(0, TT), ylab='cpue', main=sprintf("xi: %s,  zeta: %s ", xiSims[j], zetaSims[i]) )
#		#points(1:TT, hake)
#		#datGen$plotMean(add=T)
#		#if(!any(is.na(datGen$N)) & !any(datGen$N<0) & isGood(par, extra)){
#		if(!any(is.na(datGen$N)) & !any(datGen$N<0)){
#			#make data to fit
#			cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
#			#save
#			datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i]))
#			##plot
#			#points(1:TT, cpue)
#			#datGen$plotMean(add=T)
#			#datGen$plotBand()
#			
#			#
#			writeLines(sprintf("\nPID: %s, Xi: %s, Zeta: %s -> (%f, %f, %f) \n", Sys.getpid(), xiSims[j], zetaSims[i], alpha, beta, gamma))
#
#			#
#			#FIT
#			#
#		
#			#NOTE: N0Funk does not function correctly with $clone()
#        		fit = prodModel$new(
#        		        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #model
#        		        time=1:TT, catch=FtFmsy, M=M,				#constants
#        		        alpha=alpha, beta=P0, gamma=2,	#parameters
#				lalpha=log(alpha), lbeta=log(P0),	#reparameterize
#        		        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256),			#nuisance parameters
#        		        xi=xiSims[j], zeta=zetaSims[i]				#other incidentals to carry along
#        		)
#			#fit = readRDS('modsPellaFineQFixRFixP010000/fit_xi4_zeta0.35.rda')
#			fit$iterate(odeMethod)
#			#optimization
#			optAns = fit$optimize(cpue,
#			        c('lsdo', 'lalpha'), #'lq'),
#			        lower   = c(log(0.001), log(M)), #log(1e-7)),
#			        upper   = c(log(1), log(100)), #log(1e-2)),
#			        gaBoost = list(run=10, parallel=FALSE, popSize=10^3),
#				persistFor = 5,
#				fitQ    = F
#			)
#			optAns = fit$optimize(cpue,
#			        c('lsdo', 'lalpha', 'lbeta'), 
#			        lower   = c(log(0.001), log(0), log(P0/10)), 	#log(10^3)), 
#			        upper   = c(log(1), log(100), log(1.5*P0)),	#log(10^4)),
#			        gaBoost = list(run=100, parallel=FALSE, popSize=10^3),
#				persistFor = 5,
#				fitQ    = F
#			)
#			#get hessian if possible
#			tryCatch({
#				optAns = fit$optimize(cpue,
#					c('lsdo', 'lalpha', 'lbeta'), 
#                                	lower   = c(log(0.001), log(0), log(P0/10)), 
#                                	upper   = c(log(1), log(100), log(1.5*P0)), 
#					cov     = T,
#					fitQ    = F
#					#c('lsdo', 'alpha', 'beta', 'lq'),                   				        
#                                        #lower   = c(log(0.001), eps(), eps(), log(1e-7)),
#                                        #upper   = c(log(1), 10, 2, log(1e-2)),
#				)
#			}, error=function(err){
#				writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", xiSims[j], zetaSims[i]) )
#				optAns = fit$optimize(cpue,
#					c('lsdo', 'lalpha', 'lbeta'), 
#                                        lower   = c(log(0.001), log(0), log(P0/10)), 
#                                        upper   = c(log(1), log(100), log(1.5*P0)), 
#					cov     = F,
#					fitQ    = F
#					#c('lsdo', 'alpha', 'beta', 'lq'),                   				        
#                                        #lower   = c(log(0.001), eps(), eps(), log(1e-7)),
#                                        #upper   = c(log(1), 10, 2, log(1e-2)), 
#				)
#			})
#			#convenience
#			fit$alpha = exp(fit$lalpha)
#			fit$beta  = exp(fit$lbeta)
#			fit$sdo   = exp(fit$lsdo)
#			fit$q     = exp(fit$lq) 
#			#save
#			fit$save(sprintf('%s/fit_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i]))	
#			##plot
#			#fit$plotMean(add=T, col='blue')
#			#fit$plotBand(col='blue', alpha=50)
#		}
#	}
#}
##dev.off()









#
#JUNK
#

##initialize base class at a nice optimization for hake
#mod = readRDS('modAllCllhs.rds')
#datGen = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=delta, catch=catch, 
#	Cs    = mod$Cs,		#mean(catch),
#	cllhs = mod$cllhs,	#cloglog(0.5034102),
#	gamma = mod$gamma, 	#-1,
#	lq    = mod$lq,		#log(0.0004905076),
#	lsdo  = log(exp(mod$lsdo))#log(0.1160256)
#)
#datGen$model$observation = 'LN'
##
#datGen$iterate()

	##initialize base class at a nice optimization for hake
	#mod = readRDS('modAllCllhs.rds')
	#datGen = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=delta, catch=catch, 
	#	Cs    = mod$Cs,		#mean(catch),
	#	cllhs = mod$cllhs,	#cloglog(0.5034102),
	#	gamma = mod$gamma, 	#-1,
	#	lq    = mod$lq,		#log(0.0004905076),
	#	lsdo  = log(exp(mod$lsdo))#log(0.1160256)
	#)
	#datGen$model$observation = 'LN'
	##
	#datGen$iterate()




###
		##tryCatch({
                #	par = data.frame(delta=delta, cllhs=cloglog(hs), B0=B0, Bs=Bs, GA=T)
                #	CGga = ga(type    = 'real-valued',
                #	        fitness = f, pars=par,
                #	        lower   = c(0, gLow),
                #	        upper   = c(1000, 10),
		#		popSize = 100,
                #	        maxiter = 1e4,
                #	        run     = 200,
                #	        optim   = T,
                #	        monitor = F
                #	)
		#	#
                #	par$GA = F
                #	start = CGga@solution #c(300,1) #
                #	CG = multiroot(f, start, parms=par, maxiter=1e6)$root
		##}, error=function(err){ CG=c(0,0) })


##update 
		#datGen$Cs	= CG[1]	
		#datGen$gamma 	= CG[2]
		#datGen$cllhs 	= cloglog(hs)


#datGen$hs  	= cloglog(datGen$cllhs, inverse=T)



##NOTE: N0Funk does not function correctly with $clone()
#			fit = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=delta, catch=catch,
#				Cs    = CG[1],
#				gamma = -1,
#				cllhs = cloglog(hs),
#				lq    = datGen$lq,
#				lsdo  = datGen$lsdo,
#				xi    = xiSims[j],
#				zeta  = zetaSims[i]
#			)

#	#optimization
		#	optAns = fit$optimize(cpue,
		#	        c('lsdo', 'Cs', 'cllhs', 'lq'),
		#	        lower   = c(log(0.001), datGen$Cs*0.5, cloglog(0.0001), log(1e-7)),
		#	        upper   = c(log(1), datGen$Cs*2, cloglog(0.999), log(1e-2)),
		#	        gaBoost = list(run=10, parallel=FALSE, popSize=10^3)
		#	)

#	#get hessian if possible
		#	tryCatch({
		#		optAns = fit$optimize(cpue,
		#		        c('lsdo', 'Cs', 'cllhs', 'lq'),
		#		        lower   = c(log(0.001), datGen$Cs*0.5, cloglog(0.0001), log(1e-7)),
		#		        upper   = c(log(1), datGen$Cs*2, cloglog(0.999), log(1e-2)),
		#		        cov     = T
		#		)
		#	}, error=function(err){
		#		writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", xiSims[j], zetaSims[i]) )
		#		optAns = fit$optimize(cpue,
		#		        c('lsdo', 'Cs', 'cllhs', 'lq'),
		#		        lower   = c(log(0.001), datGen$Cs*0.5, cloglog(0.0001), log(1e-7)),
		#		        upper   = c(log(1), datGen$Cs*2, cloglog(0.999), log(1e-2)),
		#		        cov     = F
		#		)
		#	})


#writeLines(sprintf('\nPID: %s, Xi: %s, Zeta: %s\n', Sys.getpid(), xiSims[j], zetaSims[i]))


##
#dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
#	#linearly interpolate catches
#        ft = floor(t)
#        q  = (t-ft)
#        Cl = catch[ft]
#        Cu = catch[ft+1]
#        C = q*Cu + (1-q)*Cl
#        if(q==0){ C=Cl }
#        #
#	Fmsy = exp(lalpha)/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
#        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
#        out = R - Fmsy*C*P
#        #
#        return( list(out) )
#}
#
###
##shepSRR = function(P, alpha, beta, gamma){ alpha*P/(1+beta*P^(1/gamma)) }
#pellaSRR = function(P, alpha, beta, gamma){ alpha*P/(gamma-1)*(1-(P/beta))^(gamma-1) }
#
#
##
##P0Funk = function(alpha, beta, gamma, M){ (alpha/(M)-1)^gamma * beta^-gamma }
##exp(gamma*log(alpha/(M+ff)-1) - gamma*log(beta)) }
##PBar = function(ff, alpha, beta, gamma, M){ (alpha/(M+ff)-1)^gamma * beta^-gamma }
##PBar = function(ff, alpha, beta, gamma, M){ beta*( 1 - ((M+ff)*((gamma-1)/alpha))^(1/(gamma-1)) ) }
#PBar = function(ff, alpha, beta, gamma, M){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }
#
###
##getBeta = function(alpha, gamma, M, P0){
##        optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-P0)^2 }, method='Brent', lower=eps(), upper=10^6)$par
##}
#
##
#FMsy = function(alpha, gamma, M){
#        ##
#        #FUpper = P0 #alpha-M
#        ##root = uniroot.all(function(ff){  ((gamma-1)/alpha)^(gamma-1) - (M+ff)^(1/(gamma-1)) - (ff/(gamma-1))*(M+ff)^((2-gamma)/(gamma-1)) }, c(0, FUpper)) #1 - exp(log(ff)+log(alpha)+log(ga
#        #root = uniroot.all(function(ff){ 1 - ((M+ff)*((gamma-1)/alpha))^(1/(gamma-1)) - (ff/(gamma-1))*(M+ff)^((2-gamma)/(gamma-1))*((gamma-1)/alpha)^(1/(gamma-1)) }, c(0, FUpper))
#        ##
#        #if(length(root)<1){ return(NA) }
#        #if(length(root)>1){ return(root[1]) }
#        #
#	root = alpha/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
#        return(root)
#}
#
##
#getPar = function(xi, zeta, M){
#        #
#        gamma = 1/zeta
#        odd = (1-zeta)/zeta
#        alpha = xi*M*odd*(1-zeta)^(-odd)
#        #
#        return(c(alpha, gamma))
#}
#
###
##f = function(par, extra){
##        #unpack
##        alpha = par[1]
##        gamma = par[2]
##        #extra
##        M = extra$M
##        xi = extra$xi
##        zeta = extra$zeta
##
##        #       
##        beta = extra$P0
##
##        #par values 
##        FStar = FMsy(alpha, gamma, M)
##        PStar = PBar(FStar, alpha, beta, gamma, M)
##        #ref values
##        FSRef = xi*M
##        PSRef = zeta*extra$P0
##
##        #
##        out = c(FStar-FSRef, PStar-PSRef)
##        return( out )
##}
##
###
##howGood = function(par, extra){
##        #unpack
##        alpha = par[1]
##        gamma = par[2]
##        beta = extra$P0 #getBeta(alpha, gamma, M, extra$P0)
##        #compute
##        PZero = beta #PBar(0, alpha, beta, gamma, M)
##        fOut = f(par, extra)
##        #handel case of some numerical issue in either PBar or FMsy
##        if( length(PZero)<1 ){ return(-Inf) }
##        if( length(fOut)<2 ){ return(-Inf) }
##        #c(FStar/M-xi, PStar/PZero-zeta)
##        refComp = c(fOut[1]/extra$M, (fOut[2]+extra$P0*extra$zeta)/PZero-extra$zeta)
##        propNorm = norm(matrix(refComp, ncol=2))/norm(matrix(c(extra$xi, extra$zeta), ncol=2))
##        return( -propNorm )
##}
##isGood = function(par, extra, thresh=0.01){
##        #
##        out = (-howGood(par, extra))<thresh
##        if(is.na(out)){ out=F }
##        #
##        return( out )
##}
##
#strongRoot = function(f, par, extra, howGood, lower=c(0, 0), upper=c(2, 2), monitor=F){
#        #f      : a function computing a system of equations to find a root of
#        #par    : a vactor of starting values for computing the roots
#        #extra  : a data.frame of other incidental values required for f
#        #howGood: a function evaluating how good each root found is
#        #lower  : bounds of root search
#        #upper  : bounds of root search
#        #monitor: do you want to watch ga progress 
#        #
#        #value  : a root of f as a vector
#
#        #
#        wrap = function(x, extra){
#                #send unwanted output to /dev/null
#                capture.output(
#                #appearantly you need the <- inside capture.output else out will not be defined outside
#                out <- tryCatch({
#                        #parMR = multiroot(f, x, parms=extra, maxiter=1e6, positive=T)$root
#                        out = howGood(x, extra)
#                }, error=function(err){
#                        out = -Inf
#                }), file="/dev/null")
#                #
#                return(out)
#        }
#
#        #
#        dm = length(lower)
#        vol = prod(upper-lower)
#        gaRoot = ga(
#                type    = 'real-valued',
#                fitness = wrap,
#                extra   = extra,
#                lower   = lower,
#                upper   = upper,
#                popSize = 100*dm,
#                maxiter = 1e6,
#                run     = 200*dm,
#                optim   = T,
#                monitor = monitor#,suggestions = t(par)
#        )
#        parGA = gaRoot@solution[1,]
#        isGAGood = isGood(parGA, extra)
#        howGAGood = howGood(parGA, extra)
#        #print(isGood(parGA, extra))
#        #print(howGood(parGA, extra))
#        #print(parGA)
#        ##
#        capture.output( parMR<-multiroot(f, parGA, parms=extra, maxiter=1e6, positive=T)$root, file="/dev/null" )
#        isMRGood = isGood(parMR, extra)
#        howMRGood = howGood(parMR, extra)
#        #print( isGood(parMR, extra) )  
#        #print( howGood(parMR, extra) )
#        #print( parMR )
#
#        #
#        l = c(GA=howGAGood, MR=howMRGood)
#        who = which(l==max(l, na.rm=T))
#
#        #
#        if( isGAGood & isMRGood ){ return(list(GA=parGA, MR=parMR)[[who]]) }
#        if( isGAGood & !isMRGood ){ return(c(GA=parGA)) }
#        if( !isGAGood & isMRGood ){ return(c(MR=parMR)) }
#
#        #both are not good, but return best try anyway
#        return(list(GA=parGA, MR=parMR)[[who]])
#}

####grid for simulation
###zetaSims = seq(0.4, 0.5, 0.025) #(seq(0.1, 0.9, 0.01)) #rev(seq(0.15, 0.7, 0.01)) #rev(seq(0.1, 0.8, 0.05)) #rev(seq(0.1, 0.8, 0.01)) 	
###xiSims =   (seq(3, 4, 0.25)) #(seq(0.5, 4.5, 0.05))  #rev(seq(0.5, 3.5, 0.05)) 		#c(seq(0.5, 3.5, 0.25)) #rev(seq(0.5, 3.5, 0.05))	
##zetaSims = (seq(0.15, 0.7, 0.05))    #0.65 #seq(0.1, 0.9, 0.05) #(seq(0.1, 0.9, 0.01)) #rev(seq(0.15, 0.7, 0.01)) #rev(seq(0.1, 0.8, 0.05)) #rev(seq(0.1, 0.8, 0.01)) 	
##xiSims = rev(seq(0.5, 3.5, 0.05))  #rev(seq(0, 4, 0.5)) #(seq(0.5, 4.5, 0.05))  #rev(seq(0.5, 3.5, 0.05)) 		#c(seq(0.5, 3.5, 0.25)) #rev(seq(0.5, 3.5, 0.05))	
###zetaSims = (seq(0.45, 0.65, 0.1))    	
###xiSims = rev(seq(0.5, 4.5, 0.05))
#
##start the parameters here
#alpha = 1
#gamma = 2
#beta  = P0
#
###
##pdf('dataGrid.pdf', width=30, height=22)
##layout(matrix(1:(length(zetaSims)*length(xiSims)), nrow=length(zetaSims), ncol=length(xiSims), byrow=T))


