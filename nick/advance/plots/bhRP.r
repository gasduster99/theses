rm(list=ls())

#
library(latex2exp)

#
#
#

#
png("bhRP.png")
par(mar=c(5, 4.5, 4, 4)+0.1)
curve(1/(x+2), from=0, to=3.5, ylim=c(0.1, 0.6), lwd=3, main="BH Constrained RP-Space", xlab=TeX("F_{MSY}/M"), ylab=TeX("B_{MSY}/B_0"), cex.lab = 1.5, cex.main= 1.5)
dev.off()

#


