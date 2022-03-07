rm(list=ls())

#
library(latex2exp)

#
#
#

#
png("bhRP.png")
curve(1/(x+2), from=0, to=3.5, ylim=c(0.1, 0.6), lwd=3, main="BH Constrained RP-Space", xlab=TeX("F^*/M"), ylab=TeX("B^*/B_0"))
dev.off()

#


