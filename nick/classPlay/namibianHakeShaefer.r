rm(list=ls())

#
source('prodClass0.0.1.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, R, K, C){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #r      : growth rate given as numeric
        #K      : carrying capacity given as numeric
	#C	: catch time series
 
        #
	out = R*y*(1-y/K) - C[t]
        #
        return( list(out) )
}

#
#MAIN
#

#
cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(cpue)

#define functions
pm = prodModel$new( dNdt=dNdt, time=1:TT, N0=2709, K=2709, R=0.39, C=catch )
pm$q = 0.00045
pm$sdo 	= sd(cpue)
#pm$iterate()
#define stats model
pm$likelihood$observation = dnorm
#optimize
optOut = pm$optimize(cpue, 
	c('sdo', 'R', 'N0'), 
	lower	= c(0.01, 0, 0), 
	upper	= c(sd(cpue), 1, 1e5),  
	#gaBoost = T,
	cov	= T
)

#
#PLOT
#

#
dev.new()
plot(cpue)
lines(pm$time, pm$q*pm$N, lwd=3)
lines(pm$time, qnorm(0.025, pm$q*pm$N, pm$sdo), col='blue', lty=2)
lines(pm$time, qnorm(0.975, pm$q*pm$N, pm$sdo), col='blue', lty=2)

#
