rm(list=ls())

#
library(tgp)
library(latex2exp)

#
zetaMax = 0.6
zetaMin = 0.2
xiMax = 0.69  #3.75
xiMin = 0.1 #0.25

#
n=100

#
des = lhs(n, rbind(c(xiMin, xiMax),c(zetaMin, zetaMax)))

#
png("simplePTDesign.png")
plot(des, xlab=TeX("$F_{MSY}$"), ylab=TeX('$B_{MSY}/B_0$'), main="Latin Hypercube Design")
abline(h=0.5, lwd=3)
dev.off()



