rm(list=ls())

#
library(GA)

#
#FUNCTIONS
#

#
logline = function(q, alpha, beta, sigma, kappa, gamma){
	k = 2
	wI = 0.001
	wk = (1-kappa)*wI + kappa*wI#w_(a-1)
	#
	h  = 0
	sh = sigma*(1-h) 
        W0 = ((1-sh)*wk + (1-kappa)*sh*wI)/(1-kappa*sigma)
	rh = (wk/W0)*(1-sh)/(1-h)
	R0 = (1-(rh/alpha)^gamma)*rh/(beta*gamma) #(1-exp(gamma*(log(rh)-log(alpha))))*(rh/(beta*gamma)) #
        B0 = (W0/wk)*(R0/(1-sh)) #  (W0*Rh)/(wk*#b#10000;
        C0 = h*B0
        #
        B = matrix(NA, Y, 1)
        W = matrix(NA, Y, 1)
        R = matrix(NA, Y, 1)
        logI = matrix(NA, Y, 1)
        #
	flag = (alpha*( 1 + ( (1-kappa)*sigma )/( (1-sigma)+(1-kappa)*sigma ) * ((wI/wk) - 1) ))>(1-sigma)
	if( !flag ){ return(NA) }
	#
        R[1] = R0 #alpha*(B0-C0)*(1-beta*gamma*(B0-C0))^(1/gamma)
	B[1] = R[1] + (((1-kappa)*wI + kappa*W0)/W0) * sigma*(B0-C0) #exp(log( (1-kappa)*wI+kappa*W0 ) + log(sigma) + log(B0) - log(W0))
        W[1] = B[1]/( R[1]/wk + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1] - C[1]/2.);
        for(t in 2:Y){
		if(t>k){ R[t] = alpha*(B[t-k]-C[t-k])*(1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma)
		}else{ 	 R[t] = R0 }
       		B[t] = R[t] + (((1-kappa)*wI + kappa*W[t-1])/W[t-1]) * sigma*(B[t-1]-C[t-1]) #exp(log( (1-kappa)*
        	W[t] = B[t]/( R[t]/wk + sigma*(B[t-1]-C[t-1])/W[t-1] );
        	logI[t] = log(q) + log(B[t] - C[t]/2.);
	}
      	#print(c(sh, W0, rh, R0, B0, C0)) 
	#print(B)
	return( logI )
}

#
ll = function(par){
	#
	q = par[1]
	alpha = par[2]
	beta = par[3]
	s = par[5]
	sigma = par[4]
	#
	gamma  = -1
        kappa  = 0
	#
	logI = logline(q, alpha, beta, sigma, kappa, gamma)	
	out = sum(dnorm(log(cpue), mean=logI, sd=s, log=T))
	#
	return( out )
}

#
monitor = function(obj){
	#	
	fitness = na.exclude(obj@fitness)
	where = which(obj@fitness==max(fitness))
	writeLines(sprintf('It=%7d | Best=%3.6e | (%2.3e, %2.3e, %2.3e, %2.3e, %2.3e)', obj@iter, max(fitness), obj@population[where,1], obj@population[where,2], obj@population[where,3], obj@population[where,4], obj@population[where,5]))
	#
	plot(year, cpue)
	lines(year, exp(logline(obj@population[where,1], obj@population[where,2], obj@population[where,3], obj@population[where,4], 0, -1)))
}

#
#DATA
#

#
year = 1965:1987; 
Y = length(year)
cpue = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
C = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)

#
#OPTIMIZE
#

#
sug = rbind(
	# q, 	          alpha,     beta, 	  sigma,     s         : loglike 
	c(3.126e-04,      6.681e-02, 2.035e-01,   8.401e-01, 9.990e-1),#   
	c(exp(-7.618799), 3.002478,  0.006292,    1-0.17391, 0.3),     # 
	c(0.0004813712,   3.187019,  0.006891839, 0.844407,  0.11858)  # 1.640428e+01
	c(4.666e-04,      3.012e+00, 6.530e-03,   8.491e-01, 1.175e-01)# 1.664684e+01
	c(4.774e-04, 	  1.964e+00, 3.700e-03,   8.290e-01, 1.158e-01)# 1.694063e+01
)
nome   = c('q', 'alpha', 'beta', 'sigma', 's')
parLow = rep(.Machine$double.eps, dim(sug)[2]) #parLow[5] = 1e4
parHi  = c(1e-3, 1e1, 1e-2, 1e0, 1e0) 
##
#library(doParallel)
##workers <- rep(c("141.250.100.1", "141.250.105.3"), each = 8)
#cl = makeCluster(6, type = "PSOCK")
#registerDoParallel(cl)
#clusterExport(cl, varlist = c("Y", 'll'))
##
gaOut = ga(type="real-valued",
        fitness=ll,
        min=parLow,
        max=parHi,
        popSize=10^3,
        maxiter=10^6,
        optim=T,
	#run=10^3, 
        #parallel=cl,
	#numIslands=10,
	monitor=monitor,
	suggestions=sug,
	names=nome
)
print( gaOut@solution )
##
#plot(year, cpue) 	#q,	        alpha,    beta,     sigma,   kappa,gamma
#lines(year, exp(logline(exp(-7.618799), 3.002478, 0.006292, 1-0.17391, 0, -1)))


