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
		Blag = P0 #0
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
FF = 0.2 #0 # 0.2
k = 0.2

#
wi = 1
a0 = 2

#
alpha = 5
beta = 1
gamma = -0.59
#gamma=0 is a problem
#gamma=-1 is BH in limit

#
Rtil = alpha/(beta*(1+gamma)) * (1-gamma/(1+gamma))^(1/gamma)
N0 = Rtil/(M)#+FF)
P0 = (w(a0, wi, k)*Rtil + k*wi*N0)/(k+M)#+FF)
#P0 = 100 #10000
#N0 = 100 #100000 #P0/1000
#wi = 1 #P0/N0*10

#
TT = 45
#a0 = 2 #10 #2
dOut = dede(c(N0, P0), 1:TT, f, NULL, method="lsode") #seq(, TT, 0.1), f)
colnames(dOut) = c('time', 'N', 'B') #, 'B/N', 'Rec. Biomass', 'Ind. Growth')


#
#plot(dOut[,3]/dOut[,2], type='l')



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


