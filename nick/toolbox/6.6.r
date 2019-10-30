rm(list=ls())

#
library(RColorBrewer)

#
#FUNCTIONS
#

#return the MSY for the shaffer model
MSY = function(r, K){ r*K/4 }

#other catch
otherCatch = function(N, frac){ frac*N }

#
#GLOBAL
#

#years
Tea = 500
time = 1:(Tea+1)

#Gulf of St. Laurence
r1 = 0.15
K1 = 15234
#North Sea
r2 = 0.56
K2 = 185164

#
#PART ONE
#

#
N1msy = matrix(NA, nrow=Tea+1, ncol=1)
N2msy = matrix(NA, nrow=Tea+1, ncol=1)
#
N1other = matrix(NA, nrow=Tea+1, ncol=1)
N2other = matrix(NA, nrow=Tea+1, ncol=1)

#
N1msy[1] = K1
N2msy[1] = K2
#
N1other[1] = K1
N2other[1] = K2

#iterate
for(t in 1:Tea){
	#(1)
	N1msy[t+1] = N1msy[t] + r1*N1msy[t]*(1-N1msy[t]/K1) - MSY(r1, K1)
	N2msy[t+1] = N2msy[t] + r2*N2msy[t]*(1-N2msy[t]/K2) - MSY(r2, K2)
	
	#(2)
	N1other[t+1] = N1other[t] + r1*N1other[t]*(1-N1other[t]/K1) - otherCatch(N1other[t], 0.25)
	N2other[t+1] = N2other[t] + r2*N2other[t]*(1-N2other[t]/K2) - otherCatch(N2other[t], 0.25)
}

#
#PLOT ONE
#

#
cols = brewer.pal(9, "Set1")
#
png("./pictures/gulfOne.png")
plot(time, N1other, 'l', 
	main = "Gulf of St. Laurence (Slow Growth)",
	col  = cols[1], 
	lwd  = 3,
	ylab = "Number Individuals"
)
lines(time, N1msy, col=cols[2], lwd=3)
legend('topright', legend=c('quarter of stock', 'MSY'), lty=1, lwd=3, col=c(cols[1], cols[2]))
dev.off()
#
png("./pictures/seaOne.png")
plot(time, N2msy, 'l', 
	main = "North Sea (Fast Growth)",
	col  = cols[2], 
	lwd  = 3,
	ylab = "Number Individuals"
)
lines(time, N2other, col=cols[1], lwd=3)
legend('topright', legend=c('quarter of stock', 'MSY'), lty=1, lwd=3, col=c(cols[1], cols[2]))
dev.off()

#
#PART TWO
#

#
cRats = seq(0, 0.9, 0.1)
R = length(cRats)
#
N1maxes = matrix(NA, nrow=Tea+1, ncol=R)
N2maxes = matrix(NA, nrow=Tea+1, ncol=R)
#
r1Time = matrix(NaN, nrow=R, ncol=1)
r2Time = matrix(NaN, nrow=R, ncol=1)

#
thresh1 = 0.6*K1
thresh2 = 0.6*K2
#
C1max = 4*r1*K1/25
C2max = 4*r2*K2/25
#
N1maxes[1,] = K1/5
N2maxes[1,] = K2/5

#each c level
for(i in 1:R){
	#iterate
	for(t in 1:Tea){ 
        	#
		N1maxes[t+1,i] = N1maxes[t,i] + r1*N1maxes[t,i]*(1-N1maxes[t,i]/K1) - cRats[i]*C1max
        	N2maxes[t+1,i] = N2maxes[t,i] + r2*N2maxes[t,i]*(1-N2maxes[t,i]/K2) - cRats[i]*C2max
		#
		if(is.nan(r1Time[i]) & N1maxes[t+1,i]>thresh1){ r1Time[i]=t+1 }
		if(is.nan(r2Time[i]) & N2maxes[t+1,i]>thresh2){ r2Time[i]=t+1 }
	}
}

#
#PLOT TWO
#

png('./pictures/recovery.png')
plot(cRats*C1max, r1Time, 'l',
	main = "",
	ylab = "Recovery Time",
	xlab = "Catch Rate",
	lwd  = 3,
	xlim = c(0, max(cRats*C1max, cRats*C2max)),
	ylim = c(0, max(r1Time, r2Time)),
	col  = cols[1]
)
lines(cRats*C2max, r2Time, lwd=3, col=cols[2])
legend('bottomright', legend=c("Gulf of St. Laurence (Slow Growth)", "North Sea (Fast Growth)"), lty=1, lwd=3, col=c(cols[1], cols[2]))
dev.off()

#
#ODE PLAY
#

library(deSolve)

#
TT = 1000000
#
dNdt = function(t, y, parms){
	#
	r = parms[1]
	K = parms[2]
	C = parms[3]
	out = r*y*(1-y/K) - C*(1+sin(t))
	#
	return( list(out) )
}

#
parms = c(r1, K1, 0.9*C1max)
##
#NN = matrix(NA, nrow=TT, ncol=1)
#NN[1] = K1/5
#for(t in 1:TT){ NN[t+1]=NN[t]+dNdt(t, NN[t], parms) }
odeOut = ode(K1/5, 0:1000, dNdt, parms, "ode45")
plot(odeOut)

