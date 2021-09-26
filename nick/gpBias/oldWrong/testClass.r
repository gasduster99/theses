rm(list=ls())

#
library(VGAM)
#
source("gpClass0.0.1.r")

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

#
#MAIN
#

#
n = 30
X = matrix(c(1:n, (n+1):(2*n)), nc=2, byrow=F) #seq(-1, 10, length.out=2*n)
y = X[,1]+0.2*X[,2] #1:n #(X[,1]+0.2*X[,2])^2
fit = lm(y~X)
pp = predict(fit)
#
gp = gpModel$new( lm=fit, S2=S2, X=X, Y=y,
	v0 = 1,
	v1 = 0.1,
	s2 = 5,
	cr = -0.5
)
#
gp$fit( c('v0', 'v1', 's2', 'cr'),
	lower   = c(0, 0, 0, -1),
	upper   = c(Inf, Inf, Inf, 1),
	cov	= F
)




