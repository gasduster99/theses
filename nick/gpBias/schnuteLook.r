#rm(list=ls())

#
library(graphics)
library(latex2exp)
library(RColorBrewer)

#
SRR = function(B, alpha, beta, gamma){
        alpha*B*(1-beta*gamma*B)^(1/gamma)
}
SRR = Vectorize(SRR, 'B')

#
#MAIN
#

#
cols = brewer.pal(3, 'Set1') 

#
a = 1.5
b = 0.0001
gs = seq(-1.2-0.001,1.4, 0.2)
Cols = colorRamp(cols)( (1:length(gs))/length(gs) )
#
abit = 1000
curve(SRR(x, a, b, -1), 0, 5/b, lwd=3, ylim=c(-abit,a/b), xlim=c(0, 5/b+abit), ylab="Production", xlab="Biomass", main="Schnute Production")
i = 1
for(g in gs){
	curve(SRR(x, a, b, g), 0, 5/b, col=Cols[i], add=T)
	i = i+1
}
#text(locator(n=length(gs)-1), labels=TeX(sprintf("$\\gamma=%s$", round(gs,1)[-length(gs)])))
#dev.copy(pdf, "gammas.pdf")
#dev.off()
