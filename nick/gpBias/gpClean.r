rm(list=ls())

#
library(VGAM)
#
source("gpClass.0.0.1.r")

#
#FUNCTIONS
#

#X's are 2-vectors; L is 2x2; s2 is scalar
S2 = function(X0, X1, v0, v1, cr, s2){
        maxD = dbinorm(X0[,1], X0[,2], X0[,1], X0[,2], v0, v1, sqrt(v0*v1)*cr, log=T)
        s2*mapply(function(x01, x02, m){
                exp(dbinorm(X1[,1], X1[,2], x01, x02, v0, v1, sqrt(v0*v1)*cr, log=T)-m)
        }, X0[,1], X0[,2], maxD)#, mc.cores=1)
}



