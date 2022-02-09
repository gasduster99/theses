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

