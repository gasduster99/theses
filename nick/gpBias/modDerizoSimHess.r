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
strongRoot = function(f, par, extra, howGood, lower=c(0, 0), upper=c(2, 2), monitor=F){
        #f      : a function computing a system of equations to find a root of
        #par    : a vactor of starting values for computing the roots
        #extra  : a data.frame of other incidental values required for f
        #howGood: a function evaluating how good each root found is
        #lower  : bounds of root search
        #upper  : bounds of root search
        #monitor: do you want to watch ga progress 
        #
        #value  : a root of f as a vector

        #
        wrap = function(x, extra){
                #send unwanted output to /dev/null
                capture.output(
                #appearantly you need the <- inside capture.output else out will not be defined outside
                out <- tryCatch({
                        #parMR = multiroot(f, x, parms=extra, maxiter=1e6, positive=T)$root
                        out = howGood(x, extra)
                }, error=function(err){
                        out = -Inf
                }), file="/dev/null")
                #
                return(out)
        }

        #
        dm = length(lower)
        vol = prod(upper-lower)
        gaRoot = ga(
                type    = 'real-valued',
                fitness = wrap,
                extra   = extra,
                lower   = lower,
                upper   = upper,
                popSize = 100*dm,
                maxiter = 1e6,
                run     = 200*dm,
                optim   = T,
                monitor = monitor#,suggestions = t(par)
        )
        parGA = gaRoot@solution[1,]
        isGAGood = isGood(parGA, extra)
        howGAGood = howGood(parGA, extra)
        #print(isGood(parGA, extra))
        #print(howGood(parGA, extra))
        #print(parGA)
        ##
        capture.output( parMR<-multiroot(f, parGA, parms=extra, maxiter=1e6, positive=T)$root, file="/dev/null" )
        isMRGood = isGood(parMR, extra)
        howMRGood = howGood(parMR, extra)
        #print( isGood(parMR, extra) )  
        #print( howGood(parMR, extra) )
        #print( parMR )

        #
        l = c(GA=howGAGood, MR=howMRGood)
        who = which(l==max(l, na.rm=T))

        #
        if( isGAGood & isMRGood ){ return(list(GA=parGA, MR=parMR)[[who]]) }
        if( isGAGood & !isMRGood ){ return(parGA) }
        if( !isGAGood & isMRGood ){ return(parMR) }

        #both are not good, but return best try anyway
        return(list(GA=parGA, MR=parMR)[[who]])
}


#
dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
	#
	C = catch[t]
	#
	R = exp(lalpha)*P/(1+exp(lbeta)*P)^(gamma)
	out = R - M*P - C
	#
	return( list(out) )
}
##
#dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
#        #
#        C = catch[t]
#        #
#        R = exp(lalpha)*P/(1+exp(lbeta)*P^(1/gamma))
#        out = R - M*P - C
#        #
#        return( list(out) )
#}

#
derizoSRR = function(P, alpha, beta, gamma){ alpha*P/(1+beta*P)^(gamma) }
##
#shepSRR = function(P, alpha, beta, gamma){ alpha*P/(1+beta*P^(1/gamma)) }

#
PBar = function(ff, alpha, beta, gamma, M){ ((alpha/(M+ff))^(1/gamma)-1) * 1/beta }
##
##P0Funk = function(alpha, beta, gamma, M){ (alpha/(M)-1)^gamma * beta^-gamma }
##exp(gamma*log(alpha/(M+ff)-1) - gamma*log(beta)) }
#PBar = function(ff, alpha, beta, gamma, M){ (alpha/(M+ff)-1)^gamma * beta^-gamma }
#
getBeta = function(alpha, gamma, M, P0){
        #optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-P0)^2 }, method='Brent', lower=eps(), upper=10^6)$par
	PBar(0, alpha, 3000, gamma, 0.2)
}

#
FMsy = function(alpha, gamma, M){
        #
        root = uniroot.all(function(ff){ 1 - (alpha/(M+ff))^(-1/gamma) - (ff/gamma)*(1/(M+ff)) }, c(0, 40))
        #
	if(length(root)<1){ return(NA) }
        if(length(root)>1){ return(root[1]) }
	#
        return(root)
}
getAlpha = function(gamma, M, xi){
	optim(1, function(x){ (FMsy(x, gamma, M)-(xi*M))^2 }, method='Brent', lower=M, upper=10^3)$par
}
##
#FMsy = function(alpha, gamma, M){
#        #
#        FUpper = alpha-M
#        root = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
#        #
#        if(length(root)<1){ return(NA) }
#        if(length(root)>1){ return(root[1]) }
#        #
#        return(root)
#}

#
f = function(par, extra){
	#unpack
        alpha = par[1]
        #beta  = par[2]
        gamma = par[2]
        #extra
	M = extra$M
        xi = extra$xi
        zeta = extra$zeta
	
	#
	beta = getBeta(alpha, gamma, M, extra$P0)
	#alpha = getAlpha(gamma, M, xi)

        #par values
        FStar = FMsy(alpha, gamma, M)
        PStar = PBar(FStar, alpha, beta, gamma, M)
        #PZero = PBar(0, alpha, beta, gamma, M)
        #ref values
        FSRef = xi*M
        PSRef = zeta*extra$P0
        #PZRef = extra$P0

        #
        out = c(FStar-FSRef, PStar-PSRef) #, PZero-PZRef)
	return( out )
}
##
#f = function(par, extra){
#        #unpack
#        alpha = par[1]
#        gamma = par[2]
#        #extra
#        M = extra$M
#        xi = extra$xi
#        zeta = extra$zeta
#
#        #       
#        beta = getBeta(alpha, gamma, M, extra$P0)
#
#        #par values 
#        FStar = FMsy(alpha, gamma, M)
#        PStar = PBar(FStar, alpha, beta, gamma, M)
#        #ref values
#        FSRef = xi*M
#        PSRef = zeta*extra$P0
#
#        #
#        out = c(FStar-FSRef, PStar-PSRef)
#        return( out )
#}

#
howGood = function(par, extra){
        #unpack
        alpha = par[1]
        #beta  = par[2]
        gamma = par[2]
	
	#
	beta = getBeta(alpha, gamma, M, extra$P0)
	
        #compute
        fOut = f(par, extra)
        if( any(is.na(fOut)) ){ return(-Inf) }
	#
        PStar = fOut[2]+extra$P0*extra$zeta
        PZero = PBar(0, alpha, beta, gamma, M) #fOut[3]+extra$P0 
	#handel case of some numerical issue in either PBar or FMsy 
        if( length(fOut)<2 ){ return(-Inf) }
	if( PStar<0 ){ return(-Inf) }
        if( PZero<0 ){ return(-Inf) }	
        #c(FStar/M-xi, PStar/PZero-zeta)
        refComp = c(fOut[1]/extra$M, PStar/PZero-extra$zeta)
        propNorm = norm(matrix(refComp, ncol=2))/norm(matrix(c(extra$xi, extra$zeta), ncol=2))
        return( -propNorm )
}
isGood = function(par, extra, thresh=0.05){
        #
        out = (-howGood(par, extra))<thresh
        if(is.na(out)){ out=F }
        #
        return( out )
}
##
#howGood = function(par, extra){
#        #unpack
#        alpha = par[1]
#        gamma = par[2]
#        #compute
#        PZero = PBar(0, alpha, beta, gamma, M)
#        fOut = f(par, extra)
#        #handel case of some numerical issue in either PBar or FMsy
#        if( length(PZero)<1 ){ return(-Inf) }
#        if( length(fOut)<2 ){ return(-Inf) }
#        #c(FStar/M-xi, PStar/PZero-zeta)
#        refComp = c(fOut[1]/extra$M, (fOut[2]+extra$P0*extra$zeta)/PZero-extra$zeta)
#        propNorm = norm(matrix(refComp, ncol=2))/norm(matrix(c(extra$xi, extra$zeta), ncol=2))
#        return( -propNorm )
#}
#isGood = function(par, extra, thresh=0.01){
#        #
#        out = (-howGood(par, extra))<thresh
#        if(is.na(out)){ out=F }
#        #
#        return( out )
#}

#
#INIT
#

#
hake  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(hake)
#
M = 0.2
#
P0 = 3000

#
#SIM
#

#a place to store data
place = './modsDerizoFineQFix/'

#grid for simulation
zetaSims = rev(seq(0.15, 0.9, 0.01)) #rev(seq(0.15, 0.7, 0.01)) 	#rev(seq(0.1, 0.8, 0.05)) #rev(seq(0.1, 0.8, 0.01)) 	
xiSims =   rev(seq(0.5, 4.5, 0.05))  #rev(seq(0.5, 3.5, 0.05)) 		#c(seq(0.5, 3.5, 0.25)) #rev(seq(0.5, 3.5, 0.05))	

#start the parmaters here
gamma = 1
alpha = getAlpha(gamma, M, 2)
beta = getBeta(alpha, gamma, M, P0)

#
###start the parameters here
##alpha = 2
##gamma = 1
#beta  = getBeta(alpha, gamma, M, P0)

##
#xi = 3
#zeta = 0.35
##
#registerDoParallel(8)
#opts = list(preschedule=F)
#foreach(zeta=zetaSims, .options.multicore = opts) %dopar% {
##for(zeta in zetaSims){
#for(xi in xiSims){
#
##
#Fs = M*xi
#Ps = P0*zeta
#
##
#par = c(alpha, gamma)
#extra = data.frame(xi=xi, zeta=zeta, M=M, P0=P0)
#capture.output( MR <- multiroot(f, par, parms=extra, maxiter=1e6, positive=T), file="/dev/null" )
##
#par = MR$root
#if( !isGood(par, extra) ){
#        par = strongRoot(f, par, extra, howGood, lower=c(0, -5), upper=c(10, 10), monitor=F)
#}
#
##
#print(sprintf("XI: %s, ZETA: %s, isGood: %s", xi, zeta, isGood(par, extra)))
#
#}}

###
##pdf('dataGrid.pdf', width=30, height=22)
##layout(matrix(1:(length(zetaSims)*length(xiSims)), nrow=length(zetaSims), ncol=length(xiSims), byrow=T))
#
#
registerDoParallel(48)
opts = list(preschedule=F)
foreach(i=1:length(zetaSims), .options.multicore = opts) %dopar% {
#for(i in 1:length(zetaSims)){
	#
	datGen = prodModel$new(
	        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(0, exp(lalpha), exp(lbeta), gamma, M)}, #model
	        time=1:TT, catch=catch, M=M,                    #constants
	        alpha=alpha, beta=beta, gamma=gamma,            #parameters
		lalpha=log(alpha), lbeta=log(beta),		#reparameterize
	        lq=log(0.00049), lsdo=log(0.01160256) #log(0.1160256)		#nuisance parameters
	        #xi=xi, zeta=zeta                                #other incidentals to carry along
	)
	datGen$iterate()
	
	#
        for(j in 1:length(xiSims)){
		#
		#SKIP
		#
		
		#
		fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i])
		if( file.exists(fileFit) ){ 
			#writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
			next 
		}

		#
		#INVERT
		#
		
		#
                Fs = M*xiSims[j]
                Ps = P0*zetaSims[i]

		#
		par = c(alpha, gamma)
		extra = data.frame(xi=xiSims[j], zeta=zetaSims[i], M=M, P0=P0)
		capture.output( MR <- multiroot(f, par, parms=extra, maxiter=1e6, positive=T), file="/dev/null" )
		#
		par = MR$root
		if( !isGood(par, extra) ){
		        par = strongRoot(f, par, extra, howGood, lower=c(0, 0), upper=c(5, 15), monitor=F)
		}
		##
		#alpha = par[1]
		#beta  = par[2]
		#gamma = par[3]	
		##
		alpha = par[1]
		gamma = par[2]
		beta  = getBeta(alpha, gamma, M, P0)
		#
		#par values 
        	FStar = FMsy(alpha, gamma, M)
        	PStar = PBar(FStar, alpha, beta, gamma, M)
		PZero = PBar(0, alpha, beta, gamma, M)

		#
		writeLines(sprintf("Xi: %s, Zeta: %s -> (%f, %f, %f) -> (%f, %f) isGood: %d \n", xiSims[j], zetaSims[i], alpha, beta, gamma, FStar/M, PStar/PZero, isGood(par, extra)))

		#
		#DATA
		#
		
		#update data generator
		datGen$alpha = alpha;	datGen$lalpha = log(alpha)
		datGen$beta  = beta;	datGen$lbeta  = log(beta)
		datGen$gamma = gamma
		#convience
		datGen$sdo 	= exp(datGen$lsdo)	
		datGen$q   	= exp(datGen$lq)
		datGen$xi	= xiSims[j]
		datGen$zeta	= zetaSims[i]
		#iterate
		datGen$iterate()	
		
		##
		#datGen$printSelf()
		if( isGood(par, extra) ){ 
			#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i])) }
			datGen$printSelf()
		}

		##
		#plot(c(0,0), ylim=c(0.25, 2), xlim=c(0, TT), ylab='cpue', main=sprintf("xi: %s,  zeta: %s ", xiSims[j], zetaSims[i]) )
		#points(1:TT, hake)
		#datGen$plotMean(add=T)
		if(!any(is.na(datGen$N)) & !any(datGen$N<0) & isGood(par, extra)){
			#make data to fit
			cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
			#save
			datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i]))
			##plot
			#points(1:TT, cpue)
			#datGen$plotMean(add=T)
			#datGen$plotBand()
			
			#
			writeLines(sprintf("\nPID: %s, Xi: %s, Zeta: %s -> (%f, %f, %f) isGood: %d\n", Sys.getpid(), xiSims[j], zetaSims[i], alpha, beta, gamma, isGood(par, extra)))

			#
			#FIT
			#
		
			#NOTE: N0Funk does not function correctly with $clone()
        		fit = prodModel$new(
        		        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(0, exp(lalpha), exp(lbeta), gamma, M)}, #model
        		        time=1:TT, catch=catch, M=M,				#constants
        		        alpha=alpha, beta=getBeta(alpha, 1, M, P0), gamma=1,	#parameters
				lalpha=log(alpha), lbeta=log(getBeta(alpha, 1, M, P0)),	#reparameterize
        		        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256),			#nuisance parameters
        		        xi=xiSims[j], zeta=zetaSims[i]				#other incidentals to carry along
        		)
			#optimization
			optAns = fit$optimize(cpue,
			        c('lsdo', 'lalpha', 'lbeta'), #'lq'),
			        lower   = c(log(0.001), log(0.01), log(10^-6)), #log(1e-7)),
			        upper   = c(log(1), log(100), log(10)), #log(1e-2)),
			        gaBoost = list(run=10, parallel=FALSE, popSize=10^3)
			)
			#get hessian if possible
			tryCatch({
				optAns = fit$optimize(cpue,
					c('lsdo', 'lalpha', 'lbeta'), #'lq'),
                               	lower   = c(log(0.001), log(0.01), log(10^-6)), #log(1e-7)),
                                	upper   = c(log(1), log(100), log(10)), #log(1e-2)),
					cov     = T
					#c('lsdo', 'alpha', 'beta', 'lq'),                   				        
                                        #lower   = c(log(0.001), eps(), eps(), log(1e-7)),
                                        #upper   = c(log(1), 10, 2, log(1e-2)),
				)
			}, error=function(err){
				writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", xiSims[j], zetaSims[i]) )
				optAns = fit$optimize(cpue,
					c('lsdo', 'lalpha', 'lbeta'), #'lq'),
                                        lower   = c(log(0.001), log(0.01), log(10^-6)), #log(1e-7)),
                                        upper   = c(log(1), log(100), log(10)), #log(1e-2)),
					cov     = F
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
			fit$save(sprintf('%s/fit_xi%s_zeta%s.rda', place, xiSims[j], zetaSims[i]))	
			##plot
			#fit$plotMean(add=T, col='blue')
			#fit$plotBand(col='blue', alpha=50)
		}
	}
}
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

