rm(list=ls())

#
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
#DATA
#

#
TT = 20
M = 0.2
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)

#
zetaSims = seq(0.2, 0.6, 0.1)
xiSims = seq(0.5, 3, 0.5)

#
Fs = M*1.2 #*xiSims[1]
hs = 1-exp(-Fs)

#
B0 = 2500
Bs = B0*0.5#zetaSims[1]

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
	#
	N0  = N0Funk(delta, Cs, lhs, gamma)
	Nhs = NhsFunk(delta, Cs, lhs, gamma)
	#
	return( c(Nhs-Bs, N0-B0) )
}

#
start = c(300, -1)
par = data.frame(delta=1-exp(-M), lhs=logit(hs), B0=B0, Bs=Bs)
CG = multiroot(f, start, parms=par,
	maxiter=10^5,
	rtol=.Machine$double.eps,
	atol=.Machine$double.eps,
	ctol=.Machine$double.eps
)$root

#
datGen = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=1-exp(-M), catch=catch,
	lhs 	= logit(hs),
	Cs 	= CG[1],	
	gamma	= CG[2],
	lq	= log(0.0004),
	sdo	= 0.1
)
datGen$iterate()
datGen$plotMean()
cpue = rlnorm(TT, datGen$lq+log(datGen$N), datGen$sdo)
points(cpue)
##
#dev.new()
#plot(cpue, ylim=c(0, max(cpue)*1.1))
#datGen$plotMean(add=T, col='red')

