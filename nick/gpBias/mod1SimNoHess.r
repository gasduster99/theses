rm(list=ls())

#
library(GA)
library(boot)
library(rootSolve)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, Cs, lhs, gamma, H, delta){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #lhs    : logit of h* given as numeric
        #Cs     : C* given as numeric
        #gamma  : recruitment parameter given as numeric

        #Derivative
        C = catch[t]
        #
        hs = inv.logit(lhs)
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
N0Funk = function(delta, Cs, lhs, gamma){
        #parameterization
        hs = inv.logit(lhs)
	#
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
NhsFunk = function(delta, Cs, lhs, gamma){
        #parameterization
        hs = inv.logit(lhs)
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
	lhs   = pars$lhs
	B0    = pars$B0
	Bs    = pars$Bs
	GA    = pars$GA
	#
	N0  = N0Funk(delta, Cs, lhs, gamma)
	Nhs = NhsFunk(delta, Cs, lhs, gamma)
	
	#
	out = c(Nhs-Bs, N0-B0)
	if( GA ){ out = -log(sum(c(Nhs-Bs, N0-B0)^2)) }
	return( out )
}

#
#DATA
#

#
TT = 23
M = 0.2
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
#
datGen = readRDS('modAll.rds')

#
zetaSims = rev(seq(0.25, 0.75, 0.1))
xiSims = seq(0.1, 3.1, 0.5)+0.5

#
pdf('dataGrid.pdf', width=30, height=22)
layout(matrix(1:(length(zetaSims)*length(xiSims)), nrow=length(zetaSims), ncol=length(xiSims), byrow=T))
for(zi in 1:length(zetaSims)){
	for(xj in 1:length(xiSims)){
		#
		Fs = M*xiSims[xj]
		hs = 1-exp(-Fs)
		
		#
		B0 = 3000
		Bs = B0*zetaSims[zi]
		
		#
		writeLines(sprintf('\nxi: %s zeta: %s\n', xiSims[xj], zetaSims[zi]))	
		
		#
		par = data.frame(delta=1-exp(-M), lhs=logit(hs), B0=B0, Bs=Bs, GA=T)
		CGga = ga(type    = 'real-valued',
		        fitness = f, pars=par,
		        lower   = c(0, -10),
		        upper   = c(1000, 10),
		        #popSize = 1e3, 
		        maxiter = 1e4, 
		        run     = 200, 
		        optim   = T, 
		        #parallel= T
		        monitor = F 
		        #suggestions = start
		)
		
		#
		par$GA = F
		start = CGga@solution #c(300,1)#
		CG = multiroot(f, start, parms=par, maxiter=1e6)$root #, #rtol=.Machine$double.eps, #atol=.Machine$double.eps, #ctol=.Machine$double.eps	
		#	
		datGen$lhs   = logit(hs)
		datGen$Cs    = CG[1]
		datGen$gamma = CG[2]
		datGen$xi    = xiSims[xj]
		datGen$zeta  = zetaSims[zi] 
		#
		datGen$iterate()
		#	
		plot(c(0,0), ylim=c(0.25, 2), xlim=c(0, TT), ylab='cpue', main=sprintf("xi: %s,  zeta: %s ", xiSims[xj], zetaSims[zi]) )
		if(!any(is.na(datGen$N)) & !any(datGen$N<0)){
			#
			cpue = rlnorm(TT, datGen$lq+log(datGen$N), datGen$sdo)
			#	
			points(1:TT, cpue) #ylim=c(min(cpue)/1.2, max(cpue)*1.2), add=T)
			datGen$plotMean(add=T)
			datGen$plotBand()
			datGen$printSelf()
			datGen$save(sprintf('mods/datGen_xi%s_zeta%s.rda', xiSims[xj], zetaSims[zi]))
			
			#
			#INFERENCE
			#
			
			##
			#writeLines("\n## S ##\n")
			#mod1 = datGen$clone()
			#mod1$gamma = 1 #BH:-1, R:0, S:1
			#opt1 = mod1$optimize(cpue,
                        #        c('lq', 'sdo', 'Cs', 'lhs'),
                        #        lower   = c(log(1e-10), 0.001, datGen$Cs*0.5, logit(0.001)),
                        #        upper   = c(log(1e-2), 1, datGen$Cs*2, logit(0.5)),
                        #        #gaBoost = list('parallel'=8), #T,
                        #        cov     = F
                        #)
			#mod1$plotMean(add=T, col='green')
                        #mod1$plotBand(col='green')
			#mod1$printSelf()
			#mod1$save(sprintf('mods/shep_xi%s_zeta%s.rda', xiSims[xj], zetaSims[zi]))
			#
			##
			#writeLines("\n## R ##\n")
			#mod0 = datGen$clone()
			#mod0$gamma = .Machine$double.eps #BH:-1, R:0, S:1
			#opt0 = mod0$optimize(cpue,
                        #        c('lq', 'sdo', 'Cs', 'lhs'),
                        #        lower   = c(log(1e-10), 0.001, datGen$Cs*0.5, logit(0.001)),
                        #        upper   = c(log(1e-2), 1, datGen$Cs*2, logit(0.5)),
                        #        #gaBoost = list('parallel'=8), #T,
                        #        cov     = F
                        #)
			#mod0$plotMean(add=T, col='red')
                        #mod0$plotBand(col='red')
			#mod0$printSelf()
			#mod0$save(sprintf('mods/rick_xi%s_zeta%s.rda', xiSims[xj], zetaSims[zi]))
			
			#
			writeLines("\n## BH ##\n")
			modM1 = datGen$clone()
			modM1$gamma = -1 #BH:-1, R:0, S:1	
			optM1 = modM1$optimize(cpue,
			        c('lq', 'sdo', 'Cs', 'lhs'),
			        lower   = c(log(1e-10), 0.001, datGen$Cs*0.5, logit(0.001)),
			        upper   = c(log(1e-2), 1, datGen$Cs*2, logit(0.5)),
			        #gaBoost = list('parallel'=8), #T,
			        cov     = F
			)
			modM1$plotMean(add=T, col='blue')
			modM1$plotBand(col='blue')
			modM1$printSelf()
			modM1$save(sprintf('mods/bevH_xi%s_zeta%s.rda', xiSims[xj], zetaSims[zi]))
		}
	}
}

dev.off()







#datGen$catch = catch
#datOpt = datGen$optimize(cpue,
#        c('lq', 'sdo', 'Cs', 'lhs', 'gamma'),
#        lower   = c(log(1e-10), 0.001, datGen$Cs*0.5, logit(0.001), -5),
#        upper   = c(log(1e-2), 1, datGen$Cs*2, logit(0.5), 5),
#        gaBoost = list('parallel'=8), #T,
#        cov     = T
#)


