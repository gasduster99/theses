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
TT = 10

#define functions
pm = prodModel$new( dNdt=dNdt, time=1:TT, N0=500000, mn=0.2, mf=0.2 )
#
pm$sdo = 1
pm$likelihood$observation = dnorm
#
pm$iterate()
data = pm$N + rnorm(TT, 0, 10000)
pm$mf = 0.4
optOut = pm$optimize(data, c('sdo', 'mf'), lower=c(0, 0), upper=c(2000000, 1), cov=T)

#
#PLOT
#

#
dev.new()
plot(data)
lines(pm$time, pm$N, lwd=3)
lines(pm$time, qnorm(0.025, pm$N, pm$sdo), col='blue', lty=2)
lines(pm$time, qnorm(0.975, pm$N, pm$sdo), col='blue', lty=2)


