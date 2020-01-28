rm(list=ls())

#
source('prodClass0.0.1.r')

#
#FUNCTIONS
#

#
#dNdt
dNdt = function(t, y, alpha, beta, gamma){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #alpha  : recruitment parameter given as numeric
        #beta   : recruitment parameter given as numeric
        #gamma  : recruitment parameter given as numeric

	#Derivative
	
}

#
N0Funk = function(){
	#Virgin 
	W0 = (delta*w1 + (1-kappa)*(1-delta)*winf)/(1-kappa*(1-delta))
	r0 = (w1/W0)*delta
	R0 = (1-(r0/alpha)^gamma)*r0/(beta*gamma)
	P0 = R0/r0
	#
	return( P0 )
}

#
#
#

#
datGen = prodModel$new( dNdt=dNdt, time=1:TT )

#tuning parameter
kappa = 0
w1    = 0.001
winf  = 0.001







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
