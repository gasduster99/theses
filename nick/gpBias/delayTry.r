rm(list=ls())

#
library(deSolve)

#
#FUNCTIONS
#

#
w = function(a, wi, k){ wi*(1-exp(-k*a)) }

#
f = function(t, Y, p){ #a0, alpha, beta, gamma){
	#
	N = Y[1]
	B = Y[2]	
	#
	if( (t-a0)<1){
		Blag = 0
	}else{
		Blag = lagvalue(t-a0)[2]
	}
	#
	R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
	#
	out = c(N=NA, B=NA)
	out[1] = R - (M+FF)*N
	out[2] = w(a0, wi, k)*R + k*(wi*N-B) - (M+FF)*B
	#
	return( list(out) )
}

#
#MAIN
#

#
M = 0.2
FF = 0.2
k = 0.2

#
alpha = 5
beta = 1
gamma = -1

#
P0 = 10000
N0 = 100000 #P0/1000
wi = P0/N0*10

#
TT = 45
a0 = 2 #10 #2

#
dOut = dede(c(N0, P0), 1:TT, f, NULL) #seq(, TT, 0.1), f)

#
#JUNK
#

##
#derivs <- function(t, y, parms) {
#  if (t < 1)
#    dy <- -1
#  else
#    dy <- - lagvalue(t - 1)
#  list(c(dy))
#}
#
#dOut = dede(y=1, times=seq(0, 30, 0.1), derivs, parms=NULL)


