rm(list=ls())

#
library(boot)
#
#source('stochAgeClassR6.r')
source('ageClass0.0.1.r')

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
	FF = (Af>=t)*mf
	#FF = inv.logit(rf*(t-Af))*mf
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
am = ageModel$new( 
	dNdt=dNdt, 
	SRR=SRR, 
	time=1:TT, 
	N0=5000, 
	A=10, 
	Af=3, 
	As=3, 
	mn=0.2, 
	mf=0.2, 
	rf=1.2 
)
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
am$rDev = rep(0, TT) #exp(rnorm(TT, 0, 6)) #
am$rDev[5] = 1000#
am$sdo = 1
#
am$iterate()

#data
data = am$N + exp(rnorm(length(am$N), 0, 3))

#
am$rDev = rep(0, TT)

#optimize
optAns = am$optimize(data,
        c('rDev'),
        lower   = c(rep(0, TT)),
        upper   = c(rep(1e4, TT)),
        gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
        #method         = "Nelder-Mead", 
        cov     = T
)

##
#qlnorm(prob, log(self$q)+log(self$N), self$sdo)



##
#plot(rowSums(am$N), type='l', ylim=c(0, max(rowSums(am$N))))
#points(rowSums(data))
#matplot(am$N, type='l', add=T)
#matplot(data, add=T)
###
##am$rf = 1200
##am$iterate()


