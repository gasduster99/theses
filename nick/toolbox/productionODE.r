rm(list=ls())

library(deSolve)

#
#FUNCTIONS
#

dNdt = function(t, y, parms){
	#parms[1] : r growth rate given as a numeric
	#parms[2] : K carring capacity given as a numeric 
	
	#
	r = parms[1]
	K = parms[2]
	#
	out = r*y*(1-y/K) - r*K/4
	#
	return( list(out) )
}

#
#MAIN
#

#years
Tea = 500
time = 1:Tea

#Gulf of St. Laurence
r1 = 0.15
K1 = 15234
#North Sea
r2 = 0.56
K2 = 185164

#solve
N1 = ode(seq(8000, K1, 100), time, dNdt, c(r1, K1), "ode45")
#N2 = ode(0:K2, time, dNdt, c(r2, K2), "ode45")[,2]




