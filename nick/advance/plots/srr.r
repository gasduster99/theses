rm(list=ls())

#
library(latex2exp)
library(RColorBrewer)

#
SRR = function(B, r, K, gamma){r*B/(gamma-1)*(1-B/K)^(gamma-1)}
rGivenRandBeta = function(K, beta, R){R*(K/beta-1)/beta*(1-beta/K)^(1-K/beta)}

#
r = 1
K = 10000
Fr = 1/2
#
n=9
cols = brewer.pal(9,'Set1')
#
png("srr1.png")
gamma = 2
curve(SRR(x, r, K, gamma)/SRR(K/gamma, r, K, gamma), from=0, to=K, lwd=3, col="black", main=TeX("R(B; r, 10000, $\\gamma$)"), ylab="Production", xlab="B")
rug(K/gamma, lwd=3)
gamma = 4/3 #3/2
curve(SRR(x, r, K, gamma)/SRR(K/gamma, r, K, gamma), lwd=3, col=cols[1], add=T)
rug(K/gamma, col=cols[1], lwd=3)
gamma = 4
curve(SRR(x, r, K, gamma)/SRR(K/gamma, r, K, gamma), lwd=3, col=cols[2], add=T)
rug(K/gamma, col=cols[2], lwd=3)
dev.off()

#
png("srr1.1.png")
R = 2500
beta = K/2
curve(SRR(x, rGivenRandBeta(K, beta, R), K, K/beta), from=0, to=K, lwd=3, col="black", ylim=c(0, 1.25*R), main=TeX("P(B; r, 10000, $\\gamma$)"), ylab="Production", xlab="B")
rug(beta, lwd=3)
rs = rGivenRandBeta(K, beta, R)
gams = K/beta

beta = K/4
curve(SRR(x, rGivenRandBeta(K, beta, R), K, K/beta), lwd=3, col=cols[1], add=T)
rug(beta, col=cols[1], lwd=3)
rs = c(rs, rGivenRandBeta(K, beta, R))
gams = c(gams, K/beta)
#
beta = 3*K/4
curve(SRR(x, rGivenRandBeta(K, beta, R), K, K/beta), lwd=3, col=cols[2], add=T)
rug(beta, col=cols[2], lwd=3)
rs = c(rs, rGivenRandBeta(K, beta, R))
gams = c(gams, K/beta)
#
legend("top", 
	legend=c(
		TeX(sprintf("r=%1.1f, $\\gamma$=%1.1f", rs[1], gams[1])),
		TeX(sprintf("r=%1.1f, $\\gamma$=%1.1f", rs[2], gams[2])),
		TeX(sprintf("r=%1.1f, $\\gamma$=%1.1f", rs[3], gams[3]))
		), 
	col=c("black",cols[1:2]), 
	lwd=3, 
	lty=1
)
dev.off()

#
i = 1
n = 11
cols = brewer.pal(n,'RdYlBu')
cols[6] = 'black'
lo = 1.5
#
png("srr2.png")
curve(SRR(x, r, K, lo), from=0, to=K, lwd=3, col=cols[i], ylim=c(0, SRR(K/lo, r, K, lo)), main=TeX("P(B; 1, 10000, $\\gamma$)"), xlab="B", ylab="Production")
gams = seq(lo, 2.5, length.out=n)
for(gamma in gams){
	curve(SRR(x, r, K, gamma), lwd=3, add=T, col=rev(cols)[i])
	rug(K/gamma, col=rev(cols)[i], lwd=3)
	i = i+1
}
legend("topleft", legend=TeX(sprintf("$\\gamma$=%1.1f", gams)), col=rev(cols), lwd=3, lty=1)
dev.off()
#gamma = 2.1
#curve(SRR(x, r, K, gamma), lwd=3, col=cols[2], add=T)
#rug(K/gamma, col=cols[2], lwd=3)

#
n=9
cols = brewer.pal(9,'Set1')

#
png("srrSchaeffer.png")
curve(SRR(x, r, K, 2), from=0, to=K, lwd=3, ylab="Production", xlab="B", main="Logistic SRR and Related Quantities")
curve(Fr*x, lwd=3, col=cols[1], add=T, lty=2)
#curve(SRR(x, r, K, 2)-Fr*x, from=0, to=K, lwd=3, col=cols[2], add=T)
curve(r*x, lwd=3, col=cols[2], lty=2, add=T)
segments(K/2, SRR(K/2, r, K, 2), K/2, 0, col=cols[3], lty=3, lwd=3)
segments(K, SRR(K/2, r, K, 2), K, 0, col=cols[4], lty=3, lwd=3)
legend("topleft", legend=c("Logistic SRR", TeX("$F_{MSY}B=rB/2$"), "rB", TeX("$B_{MSY}=K/2$"), "K"), col=c("black", cols[1:4]), lwd=3, lty=c(1, 2, 2, 3, 3))
dev.off()

#
png("dBdtSchaeffer.png")
curve(SRR(x, r, K, 2)-Fr*x, from=0, to=K, lwd=3)
#curve(Fr*x, lwd=3, col=cols[1], add=T)
dev.off()

#
sm = K/2
rm = K/4
#SSRDeriso = function(B, a, b, g){ (a*B)/((1-b*g*B)^(1/g)) }
bGiven = function(sm, g){ 1/(sm*(1+g)) }
aGiven = function(sm, rm, g){ bGiven(sm, g)*rm*(1+g)^((1+g)/g) }
SSRDeriso = function(B, sm, rm, g){ (aGiven(sm, rm, g)*B)*((1-bGiven(sm, g)*g*B)^(1/g)) }
#SSRBH = function(B, sm, rm, g){ (aGiven(sm, rm, g)*B)*((1-bGiven(sm, g)*B) }
#BH
#curve(SSRDeriso(x, r, K, 0.01), from=0, to=K)
cols = brewer.pal(9,'Set1')
png("derisoSrr.png")
curve(SSRDeriso(x, sm, rm, 3), from=0, to=K*1.5, add=F, n=10000000, lwd=3, col="grey30", ylim=c(0, 3000), main=TeX("P(B; $\\theta$=$\\[\\alpha$, $1/\\beta$, $\\gamma\\]$)"), ylab="Production", xlab="B")
curve(SSRDeriso(x, sm, rm, 1), from=0, to=K, add=T, lwd=3, col=cols[1])
curve(SSRDeriso(x, sm, rm, 0.5), from=0, to=K*1.5, add=T, lwd=3, col="grey60")
curve(SSRDeriso(x, sm, rm, 0.01), from=0, to=K*1.5, add=T, lwd=3, col=cols[2])
curve(SSRDeriso(x, sm, rm, -0.5), from=0, to=K*1.5, add=T, lwd=3, col="grey85")
curve(SSRDeriso(x, sm, rm, -0.99), from=0, to=K*1.5, add=T, n=100, lwd=3, col=cols[3])
leg = c(
	TeX(sprintf("$\\theta_{B}$=( %d, %s, %s)", 50, 1/bGiven(sm, -0.99), -1)),
	TeX(sprintf("$\\theta$=(%1.1f, %1.1f, %1.1f)", aGiven(sm, rm, -0.5), 1/bGiven(sm, -0.5), -0.5)),
	TeX(sprintf("$\\theta_{R}$=(%s, %s, %s)", aGiven(sm, rm, 0), 1/bGiven(sm, 0), 0)),
	TeX(sprintf("$\\theta$=(%1.1f, %1.1f, %1.1f)", aGiven(sm, rm, 0.5), 1/bGiven(sm, 0.5), 0.5)),
	TeX(sprintf("$\\theta_{L}$=(%s, %s, %s)", aGiven(sm, rm, 1), 1/bGiven(sm, 1), 1)),
	TeX(sprintf("$\\theta$=(%1.1f, %1.1f, %1.1f)", aGiven(sm, rm, 3), 1/bGiven(sm, 3), 3))
)
legM = matrix(leg, nrow=2, ncol=3)
cc = c(cols[3], "grey85", cols[2], "grey60", cols[1], "grey30")
ccM = matrix(cc, nrow=2, ncol=3)
legend("topleft", legend=legM[1,], horiz=F, col=ccM[1,], lwd=3)
legend("topright", legend=legM[2,], horiz=F, col=ccM[2,], lwd=3)
dev.off()

#curve(SSRDeriso(x, sm, rm, -15), from=0, to=K*1.5, add=T, n=100, lwd=3, col="grey")
#curve(SSRDeriso(x, aGiven(sm, rm, 1), bGiven(sm, 1), 1), from=0, to=K, add=F)
#curve(SSRDeriso(x, aGiven(sm, rm, -0.8), bGiven(sm, -1), -1), from=0, to=K, add=T)
#curve(SSRDeriso(x, r, K, -3), from=0, to=K, add=T)
#curve(SSRDeriso(x, r, K, -1), from=0, to=K, add=T)

##BH
#curve(SSRDeriso(x, r/10000, K, 1), from=0, to=K)
##Schaefer
#curve(SSRDeriso(x, r, 1/-K, -1), from=0, to=K, add=T)



