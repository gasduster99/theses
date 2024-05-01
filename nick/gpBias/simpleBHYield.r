rm(list=ls())

library(latex2exp)

#
#FUNCTIONS
#

#
pBH = function(B, a, b){
	a*B/(1+b*B)
}

#
#MAIN
#

#
a = 1
b = 1
M = 0.2
B0 = 1/b*(a/M-1)

#
png("pBHandM.png")
f = function(x){pBH(x, a, 1)}
curve(f(x), 0, B0+1, lwd=3, ylab="Biomass", xlab="Spawning Biomass", main="BH Production and Mortality")
curve(x*M, 0, B0+1, add=T, col='red', lwd=3)
legend("bottomright", legend=c(TeX("$P_{BH}(B)$"), "MB"), col=1:2, lwd=3)
dev.off()

png("yeildBH.png")
f = function(x){pBH(x, a, 1)-M*x}
curve(f(x), 0, B0, lwd=3, ylab="Surplus Biomass", xlab="Spawning Biomass", main=TeX("$P_{BH}(B)-MB$"))
dev.off()

