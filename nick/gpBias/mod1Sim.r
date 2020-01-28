rm(list=ls())

#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
#dNdt
dNdt = function(t, y, alpha, beta, gamma, H, delta){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #alpha  : recruitment parameter given as numeric
        #beta   : recruitment parameter given as numeric
        #gamma  : recruitment parameter given as numeric

	#Derivative
	C = H[t]*y
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
TT = 20
M = 0.2
H = 1-exp(-(sin(1:TT/(pi))+1)/2)

#
datGen = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=1-exp(-M), alpha=2, beta=0.5, gamma=1, H=H, sdo=0.05 )
datGen$iterate()
#
cpue = rlnorm(TT, log(datGen$N), datGen$sdo)
#
dev.new()
plot(cpue, ylim=c(0, max(cpue)*1.1))
datGen$plotMean(add=T, col='red')

#
mod = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=1-exp(-M), alpha=1, beta=1, gamma=1, H=H, sdo=0.01 )
#optimize
optAns = mod$optimize(cpue,
        c('sdo', 'alpha', 'beta', 'gamma'),
        lower   = c(0.001, 0.001, 0.001, 0.001),
        upper   = c(0.3, 10, 5, 2),
        gaBoost = T,
        cov     = T
)
mod$plotMean(add=T)
mod$plotBand()

datGen$print()
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
