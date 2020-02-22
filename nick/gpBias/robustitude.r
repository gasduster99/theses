rm(list=ls())

#
library(geoR)
library(boot)

#
#
#

#
#
boxTrans = function(x, lambda, lambda2=0){
        if(lambda==0){return(log(x+lambda2))}
        return( ((x+lambda2)^lambda-1)/lambda )
}

#
M = 10^5
#
xi = rgamma(M, 2, 0.00005) #rlnorm(M)
zeta = 1/(xi+2)
gamma = logit(2*zeta)

#
out = boxcoxfit(gamma, lambda2=T)

dev.new()
qqnorm(gamma-mean(gamma))
qqline(gamma-mean(gamma))
dev.new()
qqnorm(boxTrans(gamma, out$lambda[1], out$lambda[2]))
qqline(boxTrans(gamma, out$lambda[1], out$lambda[2]))



