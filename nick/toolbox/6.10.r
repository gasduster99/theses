rm(list=ls())

#
library(RColorBrewer)

#
#FUNCTIONS
#

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
	if(a>=edge){ y=(fm/(fm+M))*(1-exp(-M-fm))*LtoW(vbAtoL(a, t0, k, L00), rho)*Na }
	#
	return( y )
}

#
#GLOBAL
#

#max age
A = 10
#
rho = 0.1
#
k   = 0.25
L00 = 50
t0  = 0
M   = 0.1
f  = 0.1

#
#PART ONE
#

#
N = matrix(NA, nrow=A, ncol=3)
rownames(N) = sprintf("AGE %d", 0:9)
colnames(N) = sprintf("ar%d", 2:4)
N[1,] = 500000
#
Y = matrix(0, nrow=A, ncol=3)
rownames(Y) = sprintf("AGE %d", 0:9)
colnames(Y) = sprintf("ar%d", 2:4)

#iterate
for(a in 1:(A-1)){
	#
	N[a+1,1] = knifeStep(a, N[a,1], 2, f)
	N[a+1,2] = knifeStep(a, N[a,2], 3, f)
	N[a+1,3] = knifeStep(a, N[a,3], 4, f)
	#
	Y[a+1,1] = knifeYeild(a, N[a,1], 2, f)
        Y[a+1,2] = knifeYeild(a, N[a,2], 3, f)
        Y[a+1,3] = knifeYeild(a, N[a,3], 4, f)	
}

#plot
cols = brewer.pal(9, "Set1")
png('./pictures/singleCohortN.png')
matplot(N, type='l', 
	col=cols[1:3], 
	lwd=3,
	xlab="Age",
	main="Single Cohort Dynamics"
	
)
legend('topright', legend=c("Af=2", "Af=3", "Af=4"), lty=1, lwd=3, col=c(cols[1], cols[2], cols[3]))
dev.off()
#
png('./pictures/singleCohortYeild.png')
matplot(Y, type='l', 
	col=cols[1:3], 
	lwd=3,
	xlab="Age",
	main="Single Cohort Yeild"
	
)
legend('topleft', legend=c("Af=2", "Af=3", "Af=4"), lty=1, lwd=3, col=c(cols[1], cols[2], cols[3]))
dev.off()

#
#PART TWO
#

#
FF = seq(0, 1, 0.01)
I = length(FF)
#
YY = matrix(NA, nrow=I, ncol=3)

#
for(i in 1:I){
	#
	f = FF[i]
	#
	N = matrix(NA, nrow=A, ncol=3)
	rownames(N) = sprintf("AGE %d", 0:9)
	colnames(N) = sprintf("ar%d", 2:4)
	N[1,] = 500000
	
	#
	Y = matrix(0, nrow=A, ncol=3)
	rownames(Y) = sprintf("AGE %d", 0:9)
	colnames(Y) = sprintf("ar%d", 2:4)
	
	#iterate
	for(a in 1:(A-1)){
		#
		N[a+1,1] = knifeStep(a, N[a,1], 2, f)
		N[a+1,2] = knifeStep(a, N[a,2], 3, f)
		N[a+1,3] = knifeStep(a, N[a,3], 4, f)
		#
		Y[a+1,1] = knifeYeild(a, N[a,1], 2, f)
	        Y[a+1,2] = knifeYeild(a, N[a,2], 3, f)
	        Y[a+1,3] = knifeYeild(a, N[a,3], 4, f)	
	}
	
	#
	YY[i,] = colSums(Y)
}

#
#PLOT
#

#
png('./pictures/yeildFishing.png')
matplot(FF, YY,
	type='l', 
        col=cols[1:3],
        lwd=3,
        xlab="Fishing Mortality",
        ylab="Total Yeild",
	main="Yeild with Fishing Mortality"

)
legend('topleft', legend=c("Af=2", "Af=3", "Af=4"), lty=1, lwd=3, col=c(cols[1], cols[2], cols[3]))
dev.off()


##
#A = 100
##
#Na = matrix(NA, nrow=1, ncol=3)
#colnames(N) = sprintf("ar%d", 2:4)
#Na[1,] = 500000
##
#FF = seq(0, 0.02, 0.0001)
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
#


