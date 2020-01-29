rm(list=ls())

#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
#dNdt
dNdt = function(t, y, alpha, beta, gamma, catch, delta){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #alpha  : recruitment parameter given as numeric
        #beta   : recruitment parameter given as numeric
        #gamma  : recruitment parameter given as numeric

	#Derivative
	C = catch[t]
	R = alpha*(y-C)*(1-beta*gamma*(y-C))^(1/gamma)
	out = R - delta*(y-C) - C
	#
	return( list(out) )
}

#
N0Funk = function(delta, alpha, beta, gamma){
	#Virgin 
	r0 = delta
	R0 = (1-(r0/alpha)^gamma)*r0/(beta*gamma)
	P0 = R0/r0
	#
	return( P0 )
}

#
#
#

#
cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(cpue)
#
M = 0.2

#
dev.new()
plot(cpue, ylim=c(0, max(cpue)*1.1))

#
mod = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=1-exp(-M), alpha=1.465293, beta=0.0008879473, gamma=-0.6721292, catch=catch, q=0.0004, sdo=0.1243953 )
#optimize
optAns = mod$optimize(cpue,
        c('sdo', 'alpha', 'beta', 'gamma', 'delta'),
        lower   = c(0.001, 0, -2, -2, 0.0001),
        upper   = c(0.3, 5, 2, 2, 0.999),
        gaBoost = T,
        cov     = T
)
mod$plotMean(add=T)
mod$plotBand()

mod$print()

#datGen$plotMean()
##tuning parameter
#kappa = 0
#w1    = 0.001
#winf  = 0.001







##Virgin
#R = (1-(rh/alpha)^gamma)*rh/(beta*gamma)
#
##Recruitment
#if(t>k){
#	R = alpha*(y-C[t-1])*(1-beta*gamma*(y-C[t-1]))^(1/gamma) 
#}
#out = 
#C = h*P
#out = R*y*(1-y/K) - C[t]
##
#return( list(out) )
