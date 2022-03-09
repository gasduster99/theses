rm(list=ls())

library(deSolve)

#
r=0.7
K=10000
f = function(t, B, p){list(r*B*(1-B/K))}

#
h = 1
TT = 30
C = rep(1, TT)
B = matrix(NA, nrow=TT, ncol=1)
B[1] = 1
for(i in seq(2, TT)){
	B[i] = B[i-1] + h*(r*B[i-1]*(1-B[i-1]/K)) #- K/2*C[i-1]*B[i-1])
}
odeOut = ode(B[1], 1:TT, f)#, NULL, method="euler")

png("eulerTry.png")
plot(B, xlab="Time", ylab="Biomass")
lines(odeOut)
legend("topleft", legend=c("True Dynamics", "Euler"), lty=c(1, NA), pch=c(NA, 1))
dev.off()


