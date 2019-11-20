rm(list=ls())

library(deSolve)

#
#FUNCTIONS
#

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
SRRbh = function(S, a, b){
        #S : spawning biomass
        #a : BH a parameter
        #b : BH b parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

        #
        return( a*S/(1+b*S) )	
}

#
SRRshep = function(S, a, b, c){
        #S : spawning biomass
        #a : BH a parameter
        #b : BH b parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

        #
        return( a*S/(1+(b*S)^c) )	
}

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
odeForward = function(N0, dNdt, A, TT, par, SRR, AtoL, LtoW, method="ode45"){
	#N0	: virgin biomass given as an integer
	#dNdt	: mortalitiy derivative in time given as a function
	#A	: the maximum age given as an integer
	#TT	: the time horizon given as an integer
	#par	: a list of model parameters
	#SRR	: stock-recruitment relationship given as a function
	#	
	#value 	: a matrix of the age structured population size through time.
	#		TT times given in matrix rows and A ages given in 
	#		columns.

	#
	#PARAMETERS
	#
	
	#mortality
	Af = par$Af 
	nm = par$nm 
	fm = par$fm 
	#SRR
	As    = par$As 
	bhA   = par$bhA
	bhB   = par$bhB
	shepC = par$shepC
	#growth
	rho  = par$rho  #allometric
	k    = par$k    #VB parameters
	Linf = par$Linf #
	a0   = par$a0   #

	#
        #ITERATE
        #
	
	#
	Ws = LtoW(AtoL(As:A,a0,k,Linf),rho)
	
	#initialize
        N = matrix(NA, nrow=TT, ncol=A)
        colnames(N) = sprintf("AGE %d", 1:A)
        colnames(N)[As] = sprintf("%s (As)", colnames(N)[As])
        colnames(N)[Af] = sprintf("%s (Af)", colnames(N)[Af])
        rownames(N) = sprintf("TIME %d", 1:TT)
        #solve
        N[1,] = ode(N0, time[1:A], dNdt, c(Af, nm, fm), method)[,2]
        #copy the upper diagonal over
        for(r in 2:A){ for(c in r:A){ N[r,c]=N[r-1,c] }}

	#
	i = 0
	for(t in 1:(TT-2)){
	        #Length of the diagonal which starts on row t+1 
	        D = min(A,TT-t);
	        d = ode(SRR(Ws%*%N[t,As:A], bhA, bhB, shepC), time[1:D], dNdt, c(Af, nm, fm), method)[,2]
	        for(a in 1:D){
	                N[t+a,a] = d[a]
	                i = i+1
	        }
	}
	N[TT,1] = SRR(Ws%*%N[t,As:A], bhA, bhB, shepC)
	
	#
	return( N )
}

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
eulerForward = function(N0, dNdt, A, TT, par, spawn, AtoL, LtoW, live){
        #
        #UNPACK 
        #

        #mortality
        Af = par$Af
        nm = par$nm
        fm = par$fm
        #SRR
        As  = par$As
        bhA = par$bhA
        bhB = par$bhB
	shepC = par$shepC
	#growth
        rho  = par$rho  #allometric
        k    = par$k    #VB parameters
        Linf = par$Linf #
        a0   = par$a0   #
	
        #
        #ITERATE
        #

	#
	Ws = LtoW(AtoL(As:A,a0,k,Linf),rho)

        #initialize
        N = matrix(NA, nrow=TT, ncol=A)
        colnames(N) = sprintf("AGE %d", 1:A)
        colnames(N)[As] = sprintf("%s (As)", colnames(N)[As])
        colnames(N)[Af] = sprintf("%s (Af)", colnames(N)[Af])
        rownames(N) = sprintf("TIME %d", 1:TT)
        #First row
        N[1,1] = N0
        for(a in 1:(A-1)){ 
                N[1,a+1] = live(N[1,a], a, Af, nm, fm) #N[1,a] + dNdt(a, N[1,a], c(Af, nm, fm))[[1]] #       
	}
	
        #fill the rest of the matrix
        for(t in 1:(TT-1)){
                #live and die in LA
                for(a in 1:(A-1)){
                        N[t+1,a+1] = live(N[t,a], a, Af, nm, fm) #N[t,a] + dNdt(a, N[t,a], c(Af, nm, fm))[[1]] #
                }
                #baby makin
                S = Ws %*% N[t,As:A]
                N[t+1,1] = spawn(S, bhA, bhB, shepC)
        }
	
        #
        return( N )
}

#
#GIVENS
#

#
N0 = 500000
#time
TT   = 200
time = 1:TT
#max age
A   = 10

#
par = list()
#mortality
par$Af = 3
par$nm = 0.2
par$fm = 0.2
#SRR
par$As    = 3
par$bhA   = 3
par$bhB   = 1/500
par$shepC = 0.5
#growth
par$rho  = 0.1 	#allometric
par$k    = 0.25 #VB parameters
par$Linf = 50   #
par$a0   = 0    #

#
#MAIN
#

#Node45 = odeForward(N0, dNdt, A, TT, par, SRRshep, AtoL, LtoW, method="ode45")
#Node23 = odeForward(N0, dNdt, A, TT, par, SRRshep, AtoL, LtoW, method="ode23")
#Noil = eulerForward(N0, dNdt, A, TT, par, SRRshep, AtoL, LtoW, live)

#
par$shepC = 0.5
Node0.5 = odeForward(N0, dNdt, A, TT, par, SRRshep, AtoL, LtoW, method="ode45")
#
par$shepC = 2
Node2 = odeForward(N0, dNdt, A, TT, par, SRRshep, AtoL, LtoW, method="ode45")
#
par$shepC = 1
Node1 = odeForward(N0, dNdt, A, TT, par, SRRshep, AtoL, LtoW, method="ode45")

#
#PLOT
#

library(RColorBrewer)

#
cols = brewer.pal(3, "Set1")
#
to = 1e4
png("./pictures/srrCompare.png")
curve(SRRshep(x, par$bhA, par$bhB, 1/2), 
	from=0, 
	to=to, 
	lwd=3,
	col=cols[1],
	ylab="Recruitment",
	xlab="Stock Biomass",
	main="Stock Recruitment Relationships"
)
curve(SRRshep(x, par$bhA, par$bhB, 2), 
	from=0, 
	to=to, 
	col=cols[2], 
	lwd=3,
	add=T
)
curve(SRRshep(x, par$bhA, par$bhB, 1), 
	from=0, 
	to=to, 
	col=cols[3], 
	lwd=3,
	add=T
)
legend('topleft',
        legend = c(
		sprintf("Shep(x, %s, %s, %s)", par$bhA, par$bhB, 1/2),
		sprintf("Shep(x, %s, %s, %s)", par$bhA, par$bhB, 2),
		sprintf("Shep(x, %s, %s, %s) = BH(x, %s, %s)", par$bhA, par$bhB, 1, par$bhA, par$bhB)
	),
        col    = cols,
        lwd    = 3,
        lty    = 1 
)
dev.off()

#
ageCols = colorRampPalette(brewer.pal(11, "Spectral")[c(-4,-5,-6,-7)])( 10 )
#
png("./pictures/shepAgeStructure.png")
plot(time, rowSums(Node0.5), 'l',
	#col=cols[1],
	lwd=3,
	main=sprintf("Age Structured Sheperd(x, %s, %s, %s) Model", par$bhA, par$bhB, 1/2),
	xlab="Time",
	ylab="N"
)
matplot(Node0.5,
        type = 'l',
        add  = T,
        lwd  = 1,
        col  = ageCols
)
#lines(1:A, Node0.5[1,],
#        type = 'l',
#        #col  = cols[2],
#        lwd  = 2
#)
legend('right', 
        legend = c("Total", colnames(Node0.5)),
        col    = c('black', ageCols),
        lwd    = c(3, rep(1, length(ageCols))),
        lty    = c(1, rep(1:5, 2))
)
dev.off()

##
#plot(time, rowSums(Node2), 'l',
#	col=cols[2],
#	lwd=3
#)
#matplot(Node2,
#        type = 'l',
#        add  = T,
#        lwd  = 1,
#        col  = ageCols
#)

##
#plot(time, rowSums(Node1), 'l',
#	col=cols[3],
#	lwd=3
#)
#matplot(Node1,
#        type = 'l',
#        add  = T,
#        lwd  = 1,
#        col  = ageCols
#)


#
Ws = LtoW(AtoL(par$As:A,par$a0,par$k,par$Linf),par$rho)
png("./pictures/shepWeightTime.png")
plot(time, Ws%*%t(Node0.5[,par$As:A]), 'l', 
	lwd=3,
	main="Biomass thru Time",
	ylab="Biomass",
	xlab="Time"
)
dev.off()

#
dSRRshepdS = function(S, a, b, c){
	a*( (1+(b*S)^c) - c*(b*S)^c )/(1+(b*S)^c)^2
}

#
png("./pictures/shepDerivative.png")
curve(dSRRshepdS(x, par$bhA, par$bhB, 1/2), 
	from=0, 
	to=1e4, 
	ylim=c(-0.5, 3),
	ylab="Dervative",
	xlab="S",
	main=expression(frac(a*( (1+(b*S)^c) - c*(b*S)^c ),(1+(b*S)^c)^2)),
	col=cols[1],
	lwd=3
)
curve(dSRRshepdS(x, par$bhA, par$bhB, 1), from=0, to=1e4, add=T, col=cols[3], lwd=3)
curve(dSRRshepdS(x, par$bhA, par$bhB, 2), from=0, to=1e4, add=T, col=cols[2], lwd=3)
legend('topright',
        legend = c(
                sprintf("Shep(x, %s, %s, %s)", par$bhA, par$bhB, 1/2),
                sprintf("Shep(x, %s, %s, %s)", par$bhA, par$bhB, 2),
                sprintf("Shep(x, %s, %s, %s) = BH(x, %s, %s)", par$bhA, par$bhB, 1, par$bhA, par$bhB)
        ),
        col    = cols,
        lwd    = 3,
        lty    = 1
)
dev.off()


