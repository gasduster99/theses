rm(list=ls())

#
library(tgp)
library(latex2exp)

#
#
#

#
n = 10
FF = seq(0.25, 4, length.out=n)
BB = seq(0.15, 0.7, length.out=n)
#
png("designGrid.png")
plot(-1,-1, xlim=c(min(FF), max(FF)), ylim=c(min(BB), max(BB)),
	xlab=TeX("$F_{MSY}/M$"),
	ylab=TeX("$B_{MSY}/B_0$"),
	main="LHS Elements"
)
abline(h=BB, lty=3)
abline(v=FF, lty=3)
rug(FF, lwd=3)
rug(BB, lwd=3, side=2)
text(mean(FF), BB[1], TeX("\\textit{F}"))
text(FF[1], mean(BB), TeX("\\textit{B}"))
#
i = 7
j = 6
polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), col='grey', lwd=3)
text(mean(c(FF[i], FF[i+1])), mean(c(BB[j], BB[j+1])), TeX("$\\textit{F}_i\\times\\textit{B}_j$"))
#points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
points(3, 0.5)
#
i = 1
j = 9
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 2
j = 2
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 4
j = 3
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 3
j = 1
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 5
j = 8
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 6
j = 4
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 8
j = 5
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))
#
i = 9
j = 7
#polygon(c(FF[i], FF[i+1], FF[i+1], FF[i]), c(BB[j], BB[j], BB[j+1], BB[j+1]), lwd=3)
points(runif(1, FF[i], FF[i+1]), runif(1, BB[j], BB[j+1]))

#
#rect = rbind(c(min(FF), max(FF)), c(min(BB), max(BB)))
#pp = lhs(n, rect)
#points(pp)






dev.off()
