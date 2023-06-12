#rm(list=ls())

#
library(pracma)
library(graphics)
library(latex2exp)
library(rootSolve)
library(RColorBrewer)

#
SRR = function(B, alpha, beta, gamma){
        alpha*B*(1-beta*gamma*B)^(1/gamma)
}
SRR = Vectorize(SRR, 'B')

#
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) }

#
FMsy = function(alpha, gamma, M){
        #
        if(gamma==-1){ return(sqrt(alpha*M)-M) }
        #
        FUpper = alpha-M
        root = uniroot.all(function(ff){ A(gamma, ff, M)-alpha }, c(0, FUpper))

        #root = uniroot.all(function(ff){ 1 - ff/gamma/(M+ff) - (alpha/(M+ff))^(-1/gamma) }, c(0, 
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        ##
        return(root)
}

#
msySSR = function(a, b, g, M){
	#
	o = optimize(function(B){
			SRR(B, a, b, g)-M*B
		}, 
		c(eps(), PBar(a, b, g, 0, M)),
		maximum=T
	)
	return(o$objective)
	#root = uniroot.all(function(a){ Msy(a, b, g, M)-msy }, c(M+eps(), 100))
	#return(root[1])
	#root = uniroot(function(a){ Msy(a, b, g, M)-msy }, c(M+eps(), 100))
	#return(root$root)
}

#
getAlphaSSR = function(b, g, M, msy){
	#
	f = function(a){ msySSR(a, b, g, M)-msy }
	f = Vectorize(f, 'a')
	root = uniroot.all(f, c(M+eps(), 100))
	return(root[1])
	#root = uniroot(function(a){ Msy(a, b, g, M)-msy }, c(M+eps(), 100))
	#return(root$root)
}

#
getAlpha = function(gamma, ff, M){
        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}
A = Vectorize(getAlpha, "gamma")

#
getBeta = function(alpha, gamma, M, B0){
        (1-(M/alpha)^gamma)/B0/gamma
}
B = Vectorize(getBeta, c("alpha", "gamma"))

#
#MAIN
#

#
cols = brewer.pal(4, 'Set1') 

#
a = 1.5
b = 0.0001
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

#
reg = c(Inf, 1, 0, -1, -Inf)  #rev(c(-Inf, -1, 0, 1, Inf))
gs = seq(-1.4-0.00001,1.4, 0.2)

#
b = 1
M = 0.2
#ff = 0.2
#
B0 = PBar(a, b, gs[1], 0, M)
#a = A(gs[1], ff, M)
#b = B(a, gs[1], M, B0)
png('g2.png')
curve(SRR(x*B0, a, b, gs[1]), 0, 2/b, lwd=3, ylim=c(0, a), ylab="Production", xlab="B/B0", main="Schnute Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
for(g in gs){
	#
	#ff = FMsy(a, g, M)
	B0 = PBar(a, b, g, 0, M)
	#a = A(g, ff, M)
	#b = B(a, g, M, B0)
	w = c(3, 1)[round(g, 1)%%1==0]
	curve(SRR(x*B0, a, b, g), 0, 5/b, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
	#abline(0, B0*M, col=cols[sum(g<reg)])
	#print( c(g, sum(g<reg)) )
}
dev.off()

#
png('g3.png')
curve(SRR(x*B0, a, b, gs[1])/SRR(B0, a, b, gs[1]), 0, 2/b, lwd=3, ylim=c(0, 2.5), ylab="R/R0", xlab="B/B0", main="Schnute Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
for(g in gs){
	#
	#ff = FMsy(a, g, M)
	B0 = PBar(a, b, g, 0, M)
	#a = A(g, ff, M)
	#b = B(a, g, M, B0)
	w = c(3, 1)[round(g, 1)%%1==0]
	curve(SRR(x*B0, a, b, g)/SRR(B0, a, b, g), 0, 5/b, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)
	#abline(0, B0*M, col=cols[sum(g<reg)])
	#print( c(g, sum(g<reg)) )
}
legend("topright", legend=c(TeX(sprintf("$1\\leq\\gamma$")), TeX(sprintf("$%s > \\gamma\\geq %s$",c(1, 0), c(0, -1))), TeX(sprintf("$-1>\\gamma$"))), col=cols[1:4], lwd=3)
dev.off()

#gs = seq(-1.0-0.1,1.4, 0.2)

#msy=bmsy*ff
msy = 1
#
a = getAlphaSSR(b, gs[1], M, msy)
#
png('g4.png')
curve(SRR(x*B0, a, b, gs[1])/SRR(B0, a, b, gs[1]), 0, 2/b, lwd=3, ylim=c(0, 2.5), ylab="R/R0", xlab="B/B0", main="Schnute Production", col='white', n=10000)
#curve(SRR(x*B0, a, b, -1), 0, 5/b, lwd=3, n=10000, add=T)
for(g in gs){
	#
	a = getAlphaSSR(b, g, M, msy)
	B0 = PBar(a, b, g, 0, M)
	#
	w = c(3, 1)[round(g, 1)%%1==0]
	curve(SRR(x*B0, a, b, g)/SRR(B0, a, b, g), 0, 5/b, n=10000, col=cols[sum(round(g,1)<reg)], add=T, lwd=w)	
}
dev.off()
