rm(list=ls())

#
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
        #Cs     : C* given as numeric
        #cllhs  : cloglog of h* given as numeric 
        #gamma  : recruitment parameter given as numeric
        #delta  : same as h but for natural mortality given as numeric
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
SRR = function(s, Cs, cllhs, gamma, delta){
        #s      : value os stock size
        #Cs     : C* given as numeric
        #cllhs  : cloglog of h* given as numeric 
        #gamma  : recruitment parameter given as numeric
        #delta  : same as h but for natural mortality given as numeric

        #
        hs = cloglog(cllhs, inverse=T)
        sh = (1-delta)*(1-hs)
        alpha = (1-sh)/(1-hs) * (1+(gamma*hs)/(1-sh))^(1/gamma)
        beta = (hs^2)/((1-hs)*(1-sh+gamma*hs)*Cs)
        #
        R = alpha*(s)*(1-beta*gamma*(s))^(1/gamma)
        #
        return( R )
}

#
N0Funk = function(Cs, cllhs, gamma, delta){
        #Cs     : C* given as numeric
        #cllhs  : cloglog of h* given as numeric 
        #gamma  : recruitment parameter given as numeric
        #delta  : same as h but for natural mortality given as numeric

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
        gamma = -1 #x[2]
        #
        delta = pars$delta
        cllhs = pars$cllhs
        B0    = pars$B0
        Bs    = pars$Bs
        GA    = pars$GA
        #
        N0  = N0Funk(Cs, cllhs, gamma, delta)
        #Nhs = NhsFunk(Cs, cllhs, gamma, delta)

        #
        out = N0-B0 #c(Nhs-Bs, N0-B0)
        if( GA ){ out = -log(sum(c(N0-B0)^2)) }	#-log(sum(c(Nhs-Bs, N0-B0)^2)) }
        return( out )
}

#
#DATA
#

#
M = 0.2
delta = 1-exp(-M)
#
B0 = 3000

#
#MODELS
#

#		
xi = 2		#3.4	#1.1
zeta = 0.25	#0.74	#0.6

#
bh = readRDS(sprintf('./modsFine/fit_xi%s_zeta%s.rda', xi, zeta))
sn = readRDS(sprintf('./modsFine/datGen_xi%s_zeta%s.rda', xi, zeta))

#project Schnute directly down onto BH
par = data.frame(delta=delta, cllhs=sn$cllhs, B0=B0, Bs=B0/(sn$zeta+2), GA=T)
CGga = ga(type    = 'real-valued',
        fitness = f, pars=par,
        lower   = c(0), 
        upper   = c(1000), 
        popSize = 100,
        maxiter = 1e4,
        run     = 200,
        optim   = T,
        monitor = F
)
#
par$GA = F
start = CGga@solution 
CG = multiroot(f, start, parms=par, maxiter=1e6)$root

#
bhDown = sn$clone(deep=T)
bhDown$gamma = -1
bhDown$zeta = 1/(sn$zeta+2)
bhDown$Cs = CG

#
ys = seq(0, 4.6*10^3, 0.1)
#
bhR = SRR(ys, bh$Cs, bh$cllhs, bh$gamma, bh$delta)
snR = SRR(ys, sn$Cs, sn$cllhs, sn$gamma, sn$delta)
bhDownR = SRR(ys, bhDown$Cs, bhDown$cllhs, bhDown$gamma, bhDown$delta)

#
pdf(sprintf('srrComp_xi%s_zeta%s.pdf', xi, zeta))
plot(ys, bhR, 'l', 
	lwd  = 3,
	main = sprintf('SRR Comparison for (%s, %s)', xi, zeta),
	xlab = 'Stock',
	ylab = 'Recruitment',
	ylim = c(0, max(c(bhR, snR[!is.na(snR)], bhDownR)))
)
lines(ys, snR, col='red', lwd=3)
lines(ys, bhDownR, col='blue', lwd=3)
#
legend('topleft', legend=c('B-H MLE', 'Schnute', 'B-H Down'), col=c('black', 'red', 'blue'), lwd=3)
dev.off()

