rm(list=ls())

#
library(RColorBrewer)

#
#FUNCTIONS
#

#
BHSpawn = function(S, a, b){
	#S : spawning biomass
	#a : BH a parameter
	#b : BH b parameter
	
	#
	return( a*S/(1+b*S) )
}

#allometric relationship between length and weight
LtoW = function(l, rho){ rho*l^3 }

#execute VB growth to convert ages to lengths
vbAtoL = function(a, a0, k, Linf){ Linf*(1-exp(-k*(a-a0))) }

#
knifeStep = function(a, Na, edge, fm){
	#a	: age
	#Na	: number at age
	#edge	: knife edge age of recruitment
	#fm 	: fishing mortality
	
	#assume population is not in the fishery, but if it is include fishing mortality
	Na1 = Na*exp(-M)
	if(a>=edge){ Na1 = Na*exp(-M-fm) }
	#
	return( Na1 ) 
}

#
knifeYeild = function(a, Na, edge, fm){
	#a      : age
        #Na     : number at age
        #edge   : knife edge age of recruitment
        #fm     : fishing mortality
	
	#
	y = 0
	if(a>=edge){ 
		y = (fm/(fm+M)) * (1-exp(-M-fm)) * LtoW(vbAtoL(a, t0, k, L00), rho) * Na 
	}
	#
	return( y )
}

#
#GLOBAL
#

#max time
Tea = 50
#max age
A   = 10
#allometric coefficient
rho = 0.1
#
k   = 0.25
L00 = 50
t0  = 0
M   = 0.1
f   = 0.1
#
ar  = 3 
bhA = 3
bhB = 0.002

#
N = matrix(0, nrow=Tea, ncol=A)
colnames(N) = sprintf("AGE %d", 0:9)
rownames(N) = sprintf("TIME %d", 1:Tea)
N[1,1] = 500000
#for(a in 1:(A-1)){ N[1,a+1]=knifeStep(a, N[1,a], ar, f) }

##
#Y = matrix(0, nrow=A, ncol=3)
#rownames(Y) = sprintf("AGE %d", 0:9)
#colnames(Y) = sprintf("ar%d", 2:4)

#iterate
for(t in 1:(Tea-1)){
	for(a in 1:(A-1)){
		#Live and Die in LA
		N[t+1,a+1] = knifeStep(a, N[t,a], ar, f)		
	}
	#Reproduce
	S = LtoW(vbAtoL((ar+1):A, t0, k, L00), rho)%*%N[t,(ar+1):A]
	N[t+1,1] = BHSpawn(S, bhA, bhB)
}

##plot
#cols = brewer.pal(9, "Set1")
#matplot(N, type='l', col=cols[1:3], lwd=3)
#dev.new()
#matplot(Y, type='l', col=cols[1:3], lwd=3)

##
##PART TWO
##
#
#
##
#Na = matrix(NA, nrow=1, ncol=3)
#colnames(N) = sprintf("ar%d", 2:4)
#Na[1,] = 500000
##
#FF = seq(0, 0.2, 0.01)
#FI = length(FF)
##
#Y = matrix(0, nrow=FI, ncol=3)
#rownames(Y) = sprintf("Fish Mort %f", FF)
#colnames(Y) = sprintf("ar%d", 2:4)
#
##a bunch of different fishing mortalities
#for(fi in 1:FI){
#	#iterate
#	for(a in 1:(A-1)){
#		#fish
#		Y[fi, 1] = Y[fi, 1] + knifeYeild(a, Na[1], 2, FF[fi])
#	        Y[fi, 2] = Y[fi, 2] + knifeYeild(a, Na[2], 3, FF[fi])
#	        Y[fi, 3] = Y[fi, 3] + knifeYeild(a, Na[3], 4, FF[fi])	
#		#update to the next generation
#		Na[1] = knifeStep(a, Na[1], 2, FF[fi])
#		Na[2] = knifeStep(a, Na[2], 3, FF[fi])
#		Na[3] = knifeStep(a, Na[3], 4, FF[fi])
#	}
#}
#
##plot
#dev.new()
#matplot(FF, Y, type='l', col=cols[1:3], lwd=3)



