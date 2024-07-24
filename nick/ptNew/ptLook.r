#rm(list=ls())

#
library(pracma)
library(graphics)
library(latex2exp)
library(rootSolve)
library(RColorBrewer)

#
SRR = function(P, alpha, beta, gamma){
        #alpha*B*(1-beta*gamma*B)^(1/gamma)
	alpha*P/(gamma-1)*(1-(P/beta)^(gamma-1))
}
SRR = Vectorize(SRR, 'P')

#
PBar = function(alpha, beta, gamma, ff){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }

#
FMsy = function(alpha, gamma){ alpha/gamma }

#
msySSR = function(a, b, g){
	#
	ff = FMsy(a, g)
	BStar = PBar(a, b, g, ff)
	#
	return(ff*BStar)
}

#
getAlphaSSR = function(b, g, msy){
	#
	f = function(a){ msySSR(a, b, g)-msy }
	f = Vectorize(f, 'a')
	root = uniroot.all(f, c(eps(), 100))
	return(root[1])
	#root = uniroot(function(a){ Msy(a, b, g, M)-msy }, c(M+eps(), 100))
	#return(root$root)
}

##NOTE:SCHNUTE
#getAlpha = function(gamma, ff, M){
#        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
#}
#A = Vectorize(getAlpha, "gamma")
#
##
#getBeta = function(alpha, gamma, M, B0){
#        (1-(M/alpha)^gamma)/B0/gamma
#}
#B = Vectorize(getBeta, c("alpha", "gamma"))

#
#MAIN
#

#
cols = brewer.pal(4, 'Set1') 

##
#a = 1.5
#b = 0.0001
#gs = seq(-1.2-0.001,1.4, 0.2)
##Cols = colorRamp(cols)( (1:length(gs))/length(gs) )
###
#abit = 1000
##
#curve(SRR(x, a, b, -1), 0, 5/b, lwd=3, ylim=c(-abit,a/b), xlim=c(0, 5/b+abit), ylab="Production", xlab="Biomass", main="Schnute Production")
#i = 1
#for(g in gs){
#	curve(SRR(x, a, b, g), 0, 5/b, col=Cols[i], add=T)
#	i = i+1
#}
##text(locator(n=length(gs)-1), labels=TeX(sprintf("$\\gamma=%s$", round(gs,1)[-length(gs)])))
##dev.copy(pdf, "gammas.pdf")
##dev.off()

##
#reg = c(Inf, 1, 0, -1, -Inf)  #rev(c(-Inf, -1, 0, 1, Inf))
#gs = seq(-1.4-0.00001,1.4, 0.2)

##
#b = 1
#M = 0.2
##ff = 0.2
##
#B0 = PBar(a, b, gs[1], 0, M)
##a = A(gs[1], ff, M)
##b = B(a, gs[1], M, B0)
#png('g2.png')
#curve(SRR(x*B0, a, b, gs[1]), 0, 2/b, lwd=3, ylim=c(0, a), ylab="Production", xlab="B/B0", main="Schnute Production", col='white', n=10000)
##curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
#for(g in gs){
#	#
#	#ff = FMsy(a, g, M)
#	B0 = PBar(a, b, g, 0, M)
#	#a = A(g, ff, M)
#	#b = B(a, g, M, B0)
#	w = c(3, 1)[round(g, 1)%%1==0]
#	curve(SRR(x*B0, a, b, g), 0, 5/b, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
#	#abline(0, B0*M, col=cols[sum(g<reg)])
#	#print( c(g, sum(g<reg)) )
#}
#dev.off()
#
##
#png('g3.png')
#curve(SRR(x*B0, a, b, gs[1])/SRR(B0, a, b, gs[1]), 0, 2/b, lwd=3, ylim=c(0, 2.5), ylab="R/R0", xlab="B/B0", main="Schnute Production", col='white', n=10000)
##curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
#for(g in gs){
#	#
#	#ff = FMsy(a, g, M)
#	B0 = PBar(a, b, g, 0, M)
#	#a = A(g, ff, M)
#	#b = B(a, g, M, B0)
#	w = c(3, 1)[round(g, 1)%%1==0]
#	curve(SRR(x*B0, a, b, g)/SRR(B0, a, b, g), 0, 5/b, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
#	#abline(0, B0*M, col=cols[sum(g<reg)])
#	#print( c(g, sum(g<reg)) )
#}
#legend("topright", legend=c(TeX(sprintf("$1\\leq\\gamma$")), TeX(sprintf("$%s > \\gamma\\geq %s$",c(1, 0), c(0, -1))), TeX(sprintf("$-1\\geq\\gamma$"))), col=cols[1:4], lwd=3)
#dev.off()

#gs = seq(-1.0-0.1,1.4, 0.2)
small = 0.001
gs = logseq(1+small, 5-0.2, 10) #seq(1+small, 5, 0.2) #
reg = c(-Inf, 2, Inf)

#msy=bmsy*ff
msy = 1
#
b = 10000
a = getAlphaSSR(b, gs[1], msy)
B0 = PBar(a, b, gs[1], 0)
#
png('g4PT.png')
curve(SRR(x*B0, a, b, gs[1])/SRR(B0-0.01, a, b, gs[1]), 0, 1, lwd=3, ylab="P(B)", xlab="B/B0", main="Pella-Tomlinson Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, gs[1]), 0, 1, lwd=3, ylab="P(B)", xlab="B/B0", main="Pella-Tomlinson Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
for(g in gs){
	#
	a = getAlphaSSR(b, g, msy)
	B0 = PBar(a, b, g, 0)
	print(g)
	print(a)
	#
	w = 1 #c(3, 1)[round(g, 1)%%1==0]
	curve(SRR(x*B0, a, b, g)/SRR(B0-0.01, a, b, g), 0, 1, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)	
	#curve(SRR(x*B0, a, b, g), 0, 1, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
}
g = 2
a = getAlphaSSR(b, g, msy)
B0 = PBar(a, b, g, 0)
curve(SRR(x*B0, a, b, 2)/SRR(B0-0.01, a, b, g), 0, 1, n=10000, col="black", add=T, lwd=3)
#curve(SRR(x*B0, a, b, 2), 0, 1, n=10000, col="black", add=T, lwd=3)
legend("topright", legend=c(TeX(sprintf("$\\gamma<2$")), TeX(sprintf("$\\gamma=2$")), TeX(sprintf("$\\gamma>2$"))), col=c(cols[2], "black", cols[1]), lwd=3)
dev.off()

#
png('g5PT.png')
curve(SRR(x*B0, a, b, gs[1])/SRR(B0-0.01, a, b, gs[1]), 0, 1, lwd=3, ylab="Yield", xlab="B/B0", main=TeX("$P_{pt}(B; \\theta)$"), col='white', n=10000) #ella-Tomlinson Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, gs[1]), 0, 1, lwd=3, ylab="P(B)", xlab="B/B0", main="Pella-Tomlinson Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
for(g in gs){
	#
	a = getAlphaSSR(b, g, msy)
	B0 = PBar(a, b, g, 0)
	print(g)
	print(a)
	#
	w = 1 #c(3, 1)[round(g, 1)%%1==0]
	curve(SRR(x*B0, a, b, g)/SRR(B0-0.01, a, b, g), 0, 1, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)	
	#curve(SRR(x*B0, a, b, g), 0, 1, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
}
g = 2
a = getAlphaSSR(b, g, msy)
B0 = PBar(a, b, g, 0)
curve(SRR(x*B0, a, b, 2)/SRR(B0-0.01, a, b, g), 0, 1, n=10000, col="black", add=T, lwd=3)
#curve(SRR(x*B0, a, b, 2), 0, 1, n=10000, col="black", add=T, lwd=3)
legend("topright", legend=c(TeX(sprintf("$\\gamma<2$")), TeX(sprintf("$\\gamma=2$")), TeX(sprintf("$\\gamma>2$"))), col=c(cols[2], "black", cols[1]), lwd=3)
dev.off()

#
png('g5S.png')
#
n=9
cols = brewer.pal(9,'Set1')
curve(SRR(x*B0, a, b, gs[1])/SRR(B0-0.01, a, b, gs[1]), 0, 1, lwd=3, ylab="Yield", xlab="B/B0", main=TeX("$P_{l}(B; \\theta)$"), col='white', n=10000) #ella-Tomlinson Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, gs[1]), 0, 1, lwd=3, ylab="P(B)", xlab="B/B0", main="Pella-Tomlinson Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
#for(g in gs){
#	#
#	a = getAlphaSSR(b, g, msy)
#	B0 = PBar(a, b, g, 0)
#	print(g)
#	print(a)
#	#
#	w = 1 #c(3, 1)[round(g, 1)%%1==0]
#	curve(SRR(x*B0, a, b, g)/SRR(B0-0.01, a, b, g), 0, 1, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)	
#	#curve(SRR(x*B0, a, b, g), 0, 1, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
#}
g = 2
a = getAlphaSSR(b, g, msy)
B0 = PBar(a, b, g, 0)
curve(SRR(x*B0, a, b, 2)/SRR(B0-0.01, a, b, g), 0, 1, n=10000, col="black", add=T, lwd=3)
#curve(SRR(x*B0, a, b, 2), 0, 1, n=10000, col="black", add=T, lwd=3)
#legend("topright", legend=c(TeX(sprintf("$\\gamma<2$")), TeX(sprintf("$\\gamma=2$")), TeX(sprintf("$\\gamma>2$"))), col=c(cols[2], "black", cols[1]), lwd=3)
#curve(a*x, lwd=3, col=cols[2], lty=1, add=T)
#curve(*x, lwd=3, col=cols[1], add=T, lty=1)
segments(0, 0, 1/2, SRR(B0/2, a, b, 2)/SRR(B0-0.01, a, b, g), col=cols[1], lty=1, lwd=3)
segments(1/2, SRR(B0/2, a, b, 2)/SRR(B0-0.01, a, b, g), 1/2, 0, col=cols[3], lty=1, lwd=3)
segments(1, SRR(B0/2, a, b, 2)/SRR(B0-0.01, a, b, g), 1, 0, col=cols[4], lty=1, lwd=3)
legend("topright", legend=c("Yield", TeX("$B_{0}$"), TeX("$F_{MSY}B$"), TeX("$B_{MSY}$")), col=c("black", cols[c(4, 1, 3)]), lwd=3)
dev.off()
