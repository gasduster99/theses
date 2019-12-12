rm(list=ls())

#
library(boot)
#
source('stochAgeClassR6.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, mn, mf, Af, rf){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode 
        #mn     : natural mortality given as numeric
        #mf     : fishing mortality given as numeric
	#Af     : age suseptible to fishing given as integer
	#rf	: rate of fishing selectivity

	#
	#FF = (Af>=t)*mf
	FF = inv.logit(rf*(t-Af))*mf
	out = y*(exp(-mn-FF)-1)
        #
        return( list(out) )
}

#SRR
SRR = function(S, a, b, c){
        #S : spawning biomass
        #a : BH a parameter
        #b : BH b parameter
        #c : Shepard c parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

        #
        return( a*S/(1+(b*S)^c) )
}

#
#MAIN
#

#
TT = 20

#define functions
am = ageModel$new( dNdt=dNdt, SRR=SRR, time=1:TT, N0=500000, A=10, Af=3, As=3, mn=0.2, mf=0.2, rf=1.2 )
#SRR
am$SRR_a = 3
am$SRR_b = 1/500
am$SRR_c = 2 #0.5
#AtoL
am$AtoL_Linf = 50
am$AtoL_a0   = 0
am$AtoL_k    = 0.25
#LtoW
am$LtoW_rho = 0.1
am$LtoW_psi = 3
#
am$rDev = exp(rnorm(TT, 0, 7)) #rep(0, TT)
#am$rDev[5] = 100000
#
am$iterate()

##
#am$rf = 1200
#am$iterate()


