rm(list=ls())

#
source('prodClassR6.r')

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
TT = 100

#define functions
pm = prodModel$new( dNdt=dNdt, time=1:TT, N0=500000, mn=0.2, mf=0.2 )
#generate data
pm$iterate()
data = pm$N + rnorm(TT, 0, 10000)
#initial guesses for optimization
pm$mf 	= 0.4 	#0.17	#	
pm$N0 	= 10000	#500000	#
pm$sdo 	= 1 	#7000	#
#define stats model
pm$likelihood$observation = dnorm
#optimize
optOut = pm$optimize(data, 
	c('sdo', 'mf', 'N0'), 
	lower	= c(0, 0, 0), 
	upper	= c(50000, 1, 1e6),  
	gaBoost = T,
	cov	= T
)

#
#PLOT
#

#
dev.new()
plot(data)
lines(pm$time, pm$N, lwd=3)
lines(pm$time, qnorm(0.025, pm$N, pm$sdo), col='blue', lty=2)
lines(pm$time, qnorm(0.975, pm$N, pm$sdo), col='blue', lty=2)


