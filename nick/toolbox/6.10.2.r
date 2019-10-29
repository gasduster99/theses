rm(list=ls())

#
#FUNCTIONS
#

#
live = function(Na, a, af, nm, fm){
	#Na : number at age
	#a  : age
	#af : age suseptible to fishing
	#nm : natural mortality
	#fm : fishing mortality
	#
	#value : iterate population size forward one step

	#assume population is not in the fishery, but if it is include fishing mortality
        Na1 = Na*exp(-nm)
        if(a>=af){ Na1 = Na*exp(-nm-fm) }
	#
        return( Na1 )
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
	#l 	: length
	#rho	: the allometric coefficient between length and weight
	#
	#value 	: the weight of and individual of given length

	#
	return( rho*l^3 )
}

#execute VB growth to convert ages to lengths
AtoL = function(a, a0, k, Linf){ 
	#a	: age
	#a0	: reference age
	#k	: growth constanst
	#Linf	: maximum length
	#
	#value 	: the length of an individual at age a, assuming VB growth.
		
	#
	return( Linf*(1-exp(-k*(a-a0))) ) 
}

#
#GIVENS
#

#max time
Tea = 50
#max age
A   = 10
#allometric coefficient
rho = 0.1
#VB parameters
k   = 0.25
Linf= 50
a0  = 0
nm  = 0.1
fm  = 0.5
#BH parameters
bhA = 3
bhB = 1/50000
#age suseptable to fishing
Af  = 3
#age of first spawn 
As  = 3

#
#ITERATE FORWARD
#

#spawning weights
Ws = LtoW(AtoL(As:A,a0,k,Linf),rho)

#
N = matrix(NA, nrow=Tea, ncol=A)
colnames(N) = sprintf("AGE %d", 1:10)
colnames(N)[As] = sprintf("%s (As)", colnames(N)[As])
colnames(N)[Af] = sprintf("%s (Af)", colnames(N)[Af])
rownames(N) = sprintf("TIME %d", 1:Tea)

#fill in first row
N[1,1] = 500000
for(a in 1:(A-1)){ 
	N[1,a+1] = live(N[1,a], a, Af, nm, fm) 
}

#fill the rest of the matrix
for(t in 1:(Tea-1)){
	#live and die in LA
	for(a in 1:(A-1)){
		N[t+1,a+1] = live(N[t,a], a, Af, nm, fm)
	}
	#baby makin
	S = Ws %*% N[t,As:A]
	N[t+1,1] = spawn(S, bhA, bhB)
}

#
#PLOT
#

#dev.new()
plot(rowSums(N), type='l', ylim=c(0, max(rowSums(N))), lwd=3)
matplot(N, type='l', add=T)





