rm(list=ls())

#
library(latex2exp)

#
#FUNCTIONS
#

#
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#
#MAIN
#

#
WW = 1
a0 = -1

#
#PLOT
#

#
png('vbCurves.png')

#
aS = 0.1
kappa = 10
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 10, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, WW), col='red')
segments(aS, 0, aS, ww, col='red')
#segments(0, ww, aS, ww, col='red')
points(aS, ww, pch=19, col='red')

#
aS = 1 #4 #10
kappa = 0.5 #0.2 #0.1
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 10, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, WW), col='purple', add=T)
segments(aS, 0, aS, ww, col='purple')
#segments(0, ww, aS, ww, col='blue')
points(aS, ww, pch=19, col='purple')

#
aS = 2
kappa = 0.1
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 10, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, WW), col='blue', add=T)
segments(aS, 0, aS, ww, col='blue')
#segments(0, ww, aS, ww, col='purple')
points(aS, ww, pch=19, col='blue')


legend('bottomright', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, 1, -1), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, 1, -1), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, 1, -1),2)))),
        col=c('blue', 'purple', 'red'),
        #fill=c(cols[1:3]),#, rep(NA, 3)),
        lwd=3, #c(rep(1,3), rep(2,3)), 
        #lty=c(rev(1:3)),
        cex=1
)

dev.off()

#NOTE: set width
png('vbPoster.png', width=400, height=350)

#
aS = 0.1
kappa = 10
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 10, lwd=3, xlab="Age", ylab="Biomass", main=TeX("$w(a) = w_\\infty(1-e^{-\\kappa (a+1)})$"), ylim=c(0, WW), xlim=c(-1,7), col='red')
segments(aS, 0, aS, ww, col='red')
#segments(0, ww, aS, ww, col='red')
points(aS, ww, pch=19, col='red')

#
aS = 1 #4 #10
kappa = 0.5 #0.2 #0.1
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 10, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, WW), col='purple', add=T)
segments(aS, 0, aS, ww, col='purple')
#segments(0, ww, aS, ww, col='blue')
points(aS, ww, pch=19, col='purple')

#
aS = 2
kappa = 0.1
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 10, lwd=3, xlab="Age", ylab="Biomass", ylim=c(0, WW), col='blue', add=T)
segments(aS, 0, aS, ww, col='blue')
#segments(0, ww, aS, ww, col='purple')
points(aS, ww, pch=19, col='blue')


legend('bottomright', legend=c(TeX("Slow:       $a_s=2$    $\\kappa=0.1$"), TeX("Medium: $a_s=1$    $\\kappa=0.5$"), TeX("Fast:        $a_s=0.1$ $\\kappa=10$")),
        col=c('blue', 'purple', 'red'),
        #fill=c(cols[1:3]),#, rep(NA, 3)),
        lwd=3, #c(rep(1,3), rep(2,3)), 
        #lty=c(rev(1:3)),
        cex=1
)

dev.off()





#
png('vbDrama.png')

#
aS = 10
kappa = 0.1
ww = vbGrow(aS, kappa, WW, a0)
#
curve(vbGrow(x, kappa, WW, a0), -1, 15, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, WW), col='blue')
segments(aS, 0, aS, ww, col='blue')
#segments(0, ww, aS, ww, col='blue')
points(aS, ww, pch=19, col='blue')

dev.off()

