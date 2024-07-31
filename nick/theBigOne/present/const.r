rm(list=ls())

#
library(latex2exp)
library(RColorBrewer)

#
#
#

#
cols = brewer.pal(9, 'Set1')

#
png("blankRP.png")
plot(-1, -1, 
	ylab=TeX("$B_{MSY}/B_0$"),
        xlab=TeX("$F_{MSY}/M$"),
	xlim=c(0, 3),
	ylim=c(0, 1),
	main="RP Space"
)
dev.off()

#
png("hatchRP.png")
plot(-1, -1, 
	ylab=TeX("$B_{MSY}/B_0$"),
        xlab=TeX("$F_{MSY}/M$"),
	xlim=c(0, 3),
	ylim=c(0, 1),
	main="RP Space"
)
polygon(c(0, 4, 4, 0), c(0, 0, 1, 1), angle=-45, density=10, col=cols[2], lwd=3)
dev.off()

#
png("sRP.png")
plot(-1, -1, 
        ylab=TeX("$B_{MSY}/B_0$"),
        xlab=TeX("$F_{MSY}/M$"),
        xlim=c(0, 3),
        ylim=c(0, 1),
        main="RP Space"
)
abline(h=0.5, lwd=3)
dev.off()

#
png("sbhRP.png")
plot(-1, -1, 
        ylab=TeX("$B_{MSY}/B_0$"),
        xlab=TeX("$F_{MSY}/M$"),
        xlim=c(0, 3),
        ylim=c(0, 1),
        main="RP Space"
)
abline(h=0.5, lwd=3)
curve(1/(x+2), col='red', add=T, lwd=3)
dev.off()

