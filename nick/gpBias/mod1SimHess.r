rm(list=ls())

#
library(geoR)
library(VGAM)
library(rootSolve)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, Cs, cllhs, gamma, delta, catch){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
	#Cs   	: C* given as numeric
        #cllhs 	: cloglog of h* given as numeric 
        #gamma  : recruitment parameter given as numeric
	#delta	: same as h but for natural mortality given as numeric
	#catch  : a time series of catches, given as a vector of numerics

        #Derivative
        C = catch[t]
        #
	hs = cloglog(cllhs, inverse=T)
        sh = (1-delta)*(1-hs)
        alpha = (1-sh)/(1-hs) * (1+(gamma*hs)/(1-sh))^(1/gamma)
        beta = (hs^2)/((1-hs)*(1-sh+gamma*hs)*Cs)
        #
        R = alpha*(y-C)*(1-beta*gamma*(y-C))^(1/gamma)
        out = R - delta*(y-C) - C
        #
        return( list(out) )
}

#
N0Funk = function(Cs, cllhs, gamma, delta){
	#Cs   	: C* given as numeric
        #cllhs 	: cloglog of h* given as numeric 
        #gamma  : recruitment parameter given as numeric
	#delta	: same as h but for natural mortality given as numeric

	#parameterization
	hs = cloglog(cllhs, inverse=T)
	sh = (1-delta)*(1-hs)
        alpha = (1-sh)/(1-hs) * (1+(gamma*hs)/(1-sh))^(1/gamma)
        beta = (hs^2)/((1-hs)*(1-sh+gamma*hs)*Cs)
	#Virgin 
	r0 = delta
	R0 = (1-(r0/alpha)^gamma)*r0/(beta*gamma)
	P0 = R0/r0
	#
	return( P0 )
}

#
NhsFunk = function(Cs, cllhs, gamma, delta){
        #parameterization
        hs = cloglog(cllhs, inverse=T)
        #
        sh = (1-delta)*(1-hs)
        rh = (1-sh)/(1-hs)
        #
        alpha = rh * (1+(gamma*hs)/(1-sh))^(1/gamma)
        beta  = (hs^2)/((1-hs)*(1-sh+gamma*hs)*Cs)
        #
        Rh = rh/(beta*gamma) * (1-(rh/alpha)^gamma)
        Ph = Rh/(1-sh)
        #
        return( Ph )
}

#
f = function(x, pars){
        #
        Cs    = x[1]
        gamma = x[2]
        #
        delta = pars$delta
        cllhs = pars$cllhs
        B0    = pars$B0
        Bs    = pars$Bs
        GA    = pars$GA
        #
        N0  = N0Funk(Cs, cllhs, gamma, delta)
        Nhs = NhsFunk(Cs, cllhs, gamma, delta)

        #
        out = c(Nhs-Bs, N0-B0)
        if( GA ){ out = -log(sum(c(Nhs-Bs, N0-B0)^2)) }
        return( out )
}

#
#INIT
#

#
hake  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(hake)
#
M = 0.2
delta = 1-exp(-M)

#initialize base class at a nice optimization for hake
mod = readRDS('modAllCllhs.rds')
datGen = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=delta, catch=catch, 
	Cs    = mod$Cs,		#mean(catch),
	cllhs = mod$cllhs,	#cloglog(0.5034102),
	gamma = mod$gamma, 	#-1,
	lq    = mod$lq,		#log(0.0004905076),
	lsdo  = log(exp(mod$lsdo))#log(0.1160256)
)
datGen$model$observation = 'LN'
#
datGen$iterate()

#
#SIM
#

#
zetaSims = rev(seq(0.25, 0.75, 0.1)) #0.55	#
xiSims = seq(0.5, 3.5, 0.5) #2	#

#
pdf('dataGrid.pdf', width=30, height=22)
layout(matrix(1:(length(zetaSims)*length(xiSims)), nrow=length(zetaSims), ncol=length(xiSims), byrow=T))
for(i in 1:length(zetaSims)){
        for(j in 1:length(xiSims)){
		#
		#INVERT
		#
		
		#
                Fs = M*xiSims[j]
                hs = 1-exp(-Fs)
		gLow = ((1-delta)*(1-hs)-1)/hs
                #
                B0 = 3000
                Bs = B0*zetaSims[i]

		#
                par = data.frame(delta=delta, cllhs=cloglog(hs), B0=B0, Bs=Bs, GA=T)
                CGga = ga(type    = 'real-valued',
                        fitness = f, pars=par,
                        lower   = c(0, gLow),
                        upper   = c(1000, 10),
                        maxiter = 1e4,
                        run     = 200,
                        optim   = T,
                        monitor = F
                )
		#
                par$GA = F
                start = CGga@solution #c(300,1) #
                CG = multiroot(f, start, parms=par, maxiter=1e6)$root

		#
		#DATA
		#
		
		#update 
		datGen$Cs	= CG[1]	
		datGen$gamma 	= CG[2]
		datGen$cllhs 	= cloglog(hs)
		#convience
		datGen$sdo 	= exp(datGen$lsdo)
		datGen$hs  	= cloglog(datGen$cllhs, inverse=T)
		datGen$q   	= exp(datGen$lq)
		datGen$xi	= xiSims[j]
		datGen$zeta	= zetaSims[i]
		#iterate
		datGen$iterate()
		
		#
		plot(c(0,0), ylim=c(0.25, 2), xlim=c(0, TT), ylab='cpue', main=sprintf("xi: %s,  zeta: %s ", xiSims[j], zetaSims[i]) )
		if(!any(is.na(datGen$N)) & !any(datGen$N<0)){
			#make data to fit
			cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
			#save
			datGen$save(sprintf('modsHess/datGen_xi%s_zeta%s.rda', xiSims[j], zetaSims[i]))
			#plot
			points(1:TT, cpue)
			datGen$plotMean(add=T)
			datGen$plotBand()

			#
			#FIT
			#
			
			#NOTE: N0Funk does not function correctly with $clone()
			fit = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=delta, catch=catch,
				Cs    = CG[1],
				gamma = -1,
				cllhs = cloglog(hs),
				lq    = datGen$lq,
				lsdo  = datGen$lsdo,
				xi    = xiSims[j],
				zeta  = zetaSims[i]
			)
			#optimization
			optAns = fit$optimize(cpue,
			        c('lsdo', 'Cs', 'cllhs', 'lq'),
			        lower   = c(log(0.001), datGen$Cs*0.5, cloglog(0.0001), log(1e-7)),
			        upper   = c(log(1), datGen$Cs*2, cloglog(0.999), log(1e-2)),
			        gaBoost = list(parallel=8, run=5, popSize=10^3)#, #T,
			        #cov     = T
			)
			#get hessian if possible
			tryCatch({
				optAns = fit$optimize(cpue,
				        c('lsdo', 'Cs', 'cllhs', 'lq'),
				        lower   = c(log(0.001), datGen$Cs*0.5, cloglog(0.0001), log(1e-7)),
				        upper   = c(log(1), datGen$Cs*2, cloglog(0.999), log(1e-2)),
				        #gaBoost = list(parallel=8, run=5, popSize=10^3)#, #T,
				        cov     = T
				)
			}, error=function(err){
				writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", xiSims[j], zetaSims[i]) )
				optAns = fit$optimize(cpue,
				        c('lsdo', 'Cs', 'cllhs', 'lq'),
				        lower   = c(log(0.001), datGen$Cs*0.5, cloglog(0.0001), log(1e-7)),
				        upper   = c(log(1), datGen$Cs*2, cloglog(0.999), log(1e-2)),
				        #gaBoost = list(parallel=8, run=5, popSize=10^3)#, #T,
				        cov     = F
				)
			})
			#convenience
			fit$sdo = exp(fit$lsdo)
			fit$hs  = cloglog(fit$cllhs, inverse=T)
			fit$q   = exp(fit$lq) 
			#save
			fit$save(sprintf('modsHess/fit_xi%s_zeta%s.rda', xiSims[j], zetaSims[i]))	
			#plot
			fit$plotMean(add=T, col='blue')
			fit$plotBand(col='blue', alpha=50)
		}
	}
}
dev.off()
