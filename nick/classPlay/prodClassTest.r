rm(list=ls())

#
#source('prodClassR6.r')
source('prodClass0.0.1.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, mn, mf){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #mn     : natural mortality given as numeric
        #mf     : fishing mortality given as numeric
 
        #assume population is not in the fishery, but if it is include fishing mortality 
	out = y*(exp(-mn-mf)-1) 
        #
        return( list(out) )
}

#
#MAIN
#

#
data = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
TT = length(data)

#define functions
pm = prodModel$new( dNdt=dNdt, time=1:TT, N0=500000, mn=0.2, mf=0.2 )
#generate data
pm$iterate()
rsd = 10000
#data = pm$N + rnorm(TT, 0, rsd)
pm$sdo = rsd
#define stats model
pm$likelihood$observation = dnorm
##optimize
#optOut = pm$optimize(data, 
#	c('sdo', 'mf', 'N0'), 
#	lower	= c(0, 0, 0), 
#	upper	= c(sd(data), 1, 1e6),  
#	gaBoost = list(popSize=1e3, maxiter=1e4, run=20),
#	cov	= T
#)
#
##
##PLOT
##
#
##
#dev.new()
#plot(data)
#lines(pm$time, pm$N, lwd=3)
#lines(pm$time, qnorm(0.025, pm$N, pm$sdo), col='blue', lty=2)
#lines(pm$time, qnorm(0.975, pm$N, pm$sdo), col='blue', lty=2)


