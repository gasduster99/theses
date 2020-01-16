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
        #r      : growth rate given as numeric
        #K/N0   : carrying capacity given as numeric
        #C      : catch time series

	#Virgin
	R = (1-(rh/alpha)^gamma)*rh/(beta*gamma)
	
        #Recruitment
	if(t>k){
		R = alpha*(y-C[t-1])*(1-beta*gamma*(y-C[t-1]))^(1/gamma) 
	}
	out = 
	C = h*P
        out = R*y*(1-y/K) - C[t]
        #
        return( list(out) )
}

#
#
#

datGen = prodModel$new( dNdt=dNdt, time=1:TT, N0=3955.556, K=345.843, R=0.65, C=catch )

