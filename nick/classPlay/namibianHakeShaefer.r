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
        #K/N0   : carrying capacity given as numeric
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

#
pmLN = prodModel$new( dNdt=dNdt, time=1:TT, N0=3955.556, K=2345.843, R=0.465, C=catch )
pmLN$q = 0.00045
pmLN$sdo  = 0.1 
pmLN$model$observation = 'LN'
#optimize
optAns = pmLN$optimize(cpue, 
	c('sdo', 'R', 'K'), 
	lower	= c(0.001, 0, 0), 
	upper	= c(0.3, 1, 1e5),  
	gaBoost = F, 
	#list(maxiter=1e3, run=50, popSize=1e5), #T, #
	#method 	= "Nelder-Mead", 
	cov	= T	
)

#
dev.new()
plot(cpue)
pmLN$plot(add=T)

#
#
pmN = prodModel$new( dNdt=dNdt, time=1:TT, N0=3955.556, K=3955.556, R=0.38, C=catch )
pmN$q = 0.00045
pmN$sdo = 0.1
pmN$model$observation = 'N'
#
optAns = pmN$optimize(cpue, 
	c('sdo', 'R', 'K'), 
	lower	= c(0.001, 0, 0), 
	upper	= c(0.3, 1, 1e5),  
	gaBoost = T, 
	cov	= T	
)

#
pmN$plot(add=T, col='red')
