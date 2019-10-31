rm(list=ls())

library(deSolve)

#
#FUNCTIONS
#

#
dNdt = function(t, y, params){ 
        #params[1] : af age suseptible to fishing
        #params[2] : nm natural mortality
        #params[3] : fm fishing mortality

	#unpack
	af = params[1]
	nm = params[2]
	fm = params[3]		
        #assume population is not in the fishery, but if it is include fishing mortality
        out = y*(exp(-nm)-1)
        if(t>=af){ out = y*(exp(-nm-fm)-1) }
        #
        return( list(out) )
}

#
spawn = function(S, a, b){
        #S : spawning biomass
        #a : BH a parameter
        #b : BH b parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

        #
        return( a*S/(1+b*S) )
}

#allometric relationship between length and weight
LtoW = function(l, rho){ 
        #l      : length
        #rho    : the allometric coefficient between length and weight
        #
        #value  : the weight of and individual of given length

        #
        return( rho*l^3 )
}

#execute VB growth to convert ages to lengths
AtoL = function(a, a0, k, Linf){ 
        #a      : age
        #a0     : reference age
        #k      : growth constanst
        #Linf   : maximum length
        #
        #value  : the length of an individual at age a, assuming VB growth.

        #
        return( Linf*(1-exp(-k*(a-a0))) )
}




#
#GIVENS
#

#time
Tea  = 2000
time = 1:Tea
#max age
A   = 10
#age suseptable to fishing
Af  = 3
#age of first spawn 
As  = 3
#allometric coefficient
rho = 0.1
#VB parameters
k   = 0.25
Linf= 50
a0  = 0
nm  = 0.2
fm  = 0.4
#BH parameters
bhA = 3
bhB = 1/500

#
#ITERATE
#

#spawning weights
Ws = LtoW(AtoL(As:A,a0,k,Linf),rho)

#
N = matrix(NA, nrow=Tea, ncol=A)
colnames(N) = sprintf("AGE %d", 1:10)
colnames(N)[As] = sprintf("%s (As)", colnames(N)[As])
colnames(N)[Af] = sprintf("%s (Af)", colnames(N)[Af])
rownames(N) = sprintf("TIME %d", 1:Tea)


#solve
N[1,] = ode(500000, time[1:A], dNdt, c(Af, nm, fm), "ode45")[,2]
#copy the upper diagonal over
for(r in 2:A){ for(c in r:A){ N[r,c]=N[r-1,c] }}

#
for(t in 1:(Tea-2)){
	#Length of the diagonal which starts on row t+1 
	D = min(A,Tea-t) 
	d = ode(spawn(Ws%*%N[t,As:A], bhA, bhB), time[1:D], dNdt, c(Af, nm, fm), "ode45")[,2]
	for(a in 1:D){
		N[t+a,a] = d[a]
	} 
}
N[Tea,1] = spawn(Ws%*%N[t,As:A], bhA, bhB)


##
#live = function(Na, a, af, nm, fm){
#        #Na : number at age
#        #a  : age
#        #af : age suseptible to fishing
#        #nm : natural mortality
#        #fm : fishing mortality
#        #
#        #value : iterate population size forward one step
#
#        #assume population is not in the fishery, but if it is include fishing mortality
#        Na1 = Na*exp(-nm)
#        if(a>=af){ Na1 = Na*exp(-nm-fm) }
#        #
#        return( Na1 )
#}
#
##fill in first row
#N[1,1] = 500000
#for(a in 1:(A-1)){ 
#        N[1,a+1] = live(N[1,a], a, Af, nm, fm)
#}
#
##fill the rest of the matrix
#for(t in 1:(Tea-1)){
#        #live and die in LA
#        for(a in 1:(A-1)){
#                N[t+1,a+1] = live(N[t,a], a, Af, nm, fm)
#                #Y[t, a] = yeild(a, N[t,a], LtoW(AtoL(a, a0, k, Linf), rho), Af, nm, fm)
#        }
#        #Y[t, A] = yeild(A, N[t,A], LtoW(AtoL(A, a0, k, Linf), rho), Af, nm, fm)
#        #baby makin
#        S = Ws %*% N[t,As:A]
#        N[t+1,1] = spawn(S, bhA, bhB)
#}


