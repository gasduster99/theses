rm(list=ls())

#
library(geoR)
library(boot)

#
#
#

#
boxTrans = function(x, lambda, lambda2){
	if(lambda==0){return(log(x))}
	return( (x^lambda-1)/lambda )
}

#
M = 10^5
#
xi = rexp(M, 0.1)#rgamma(M, 1, 0.00005) #rlnorm(M)
zeta = 1/(xi+2)
gamma = logit(2*zeta)

#
ins = gamma-min(gamma)+.Machine$double.eps
out = boxcoxfit(ins) #, lambda2=T)

dev.new()
qqnorm(ins)
dev.new()
qqnorm(boxTrans(ins, out$lambda))




