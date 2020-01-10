rm(list=ls())

library(GA)

#
#FUNCTIONS
#

#
wk = function(lamb, kapp){
	wkk = 0;
        for(a in 1:k){
                wkk = lamb + kapp*wkk;
        }
	return(wkk)
}  

#
ll = function(par, x, ctch){ 
	#
	sigma  = par[3]
	lambda = par[4]
	kappa  = par[5]
	alpha  = par[6]
	beta   = par[7]
	gamma  = par[8]
	s      = par[1]
	q      = par[2]
	#
	end = 23
	#
	B0 = 2700
	R0 = 1
	C0 = 0
	W0 = 2700
	#
	B = matrix(NA, nrow=end, ncol=1)
	W = matrix(NA, nrow=end, ncol=1)
	logR = matrix(NA, nrow=end, ncol=1)
	logI = matrix(NA, nrow=end, ncol=1)	
	#
	B[1] = R0 + sigma * (lambda+kappa*W0) * (B0-C0) / W0;
        logR[1] = log(alpha) + log(B0-C0) + log(1-beta*gamma*(B0-C0))/gamma;
        W[1] = B[1]/( exp(logR[1])/wk(lambda, kappa) + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1]-ctch[1]/2.);
	for(t in 2:end){
		B[t] = exp(logR[t-1]) + sigma * (lambda+kappa*W[t-1]) * (B[t-1]-ctch[t-1]) / W[t-1];
        	logR[t] = log(alpha) + log(B[t-1]-ctch[t-1]) + log(1-beta*gamma*(B[t-1]-ctch[t-1]))/gamma;
        	W[t] = B[t]/( exp(logR[t])/wk(lambda, kappa) + sigma*(B[t-1]-ctch[t-1])/W[t-1] );
        	logI[t] = log(q) + log(B[t]-ctch[t]/2.);
	}
	#
        out = sum(dnorm(log(x), mean=logI, sd=s, log=T))
        #
        return(out)
}

#
#MAIN
#

#data
year = 1965:1987; Y = length(year)
y = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
C = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
k = 1
#optimization
s0   = 0.1
q0   = 0.001
sig0 = 0.5
lam0 = 1
kap0 = 1
alp0 = 1
bet0 = 0.5
gam0 = 0.5
par    = c(s0, q0, sig0, lam0, kap0, alp0, bet0, gam0)
parLow = rep(.Machine$double.eps, length(par))
parHi  = c(1e0, 1e0, 1e0, 1e4, 1e4, 1e0, 1e0, 1e4)
gaOut  = ga(type="real-valued",
        fitness=ll, x=y, ctch=C,
        min=parLow,
        max=parHi,
        popSize=10^3,
        maxiter=10^3,
        optim=T#,
        #parallel=T
)
colnames(gaOut@solution) = c('s', 'q', 'sigma', 'lambda', 'kappa', 'alpha', 'beta', 'gamma')
print( gaOut@solution )
#
save.image('gaOut.RData')
