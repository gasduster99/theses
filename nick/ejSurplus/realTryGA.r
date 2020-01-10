rm(list=ls())

#
library(GA)

#
#FUNCTIONS
#

#
ll = function(par){
	#
	q = par[1]
	alpha = par[2]
	#beta = par[3]
	s = par[3]
	sigma = par[4]
	w = par[5]
	b = par[6]
	#lambda = par[5]
	#kappa = par[6]
	#omega = par[7]
	#
	gamma  = -1;
        #sigma  = 1;
        lambda = 0.000731 #1;
        kappa  = 1.018862 #1;
        omega  = 0.0007;
        k = 1;
        R0 = 1;
        W0 = w#1/10000;
        B0 = b#10000;
        C0 = 0;
	#
	B = matrix(NA, Y, 1)
	W = matrix(NA, Y, 1)
	R = matrix(NA, Y, 1)
	logR = matrix(NA, Y, 1)
	logI = matrix(NA, Y, 1)
	#//logR[1] = log(alpha) + log(B0-C0) + log(1-beta*gamma*(B0-C0))/gamma; //s 
	R[1] = alpha*(B0-C0) #//exp(logR[1]);
        B[1] = R[1] + sigma * (lambda+kappa*W0) * (B0-C0) / W0;
        W[1] = B[1]/( R[1]/omega + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1] - C[1]/2.);
        for(t in 2:Y){
		R[t] = alpha*(B[t-k]-C[t-k]) #* (1-beta*(B[t-k]-C[t-k]))#* (1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma) #exp(logR[t]);
                B[t] = R[t] + sigma * (lambda+kappa*W[t-1]) * (B[t-1]-C[t-1]) / W[t-1];
                W[t] = B[t]/( R[t]/omega + sigma*(B[t-1]-C[t-1])/W[t-1] );
                logI[t] = log(q) + log(B[t] - C[t]/2.);
        }
	#print(logI)
	#print(B)
	#print(R)
	#print(W)	
	out = sum(dnorm(log(cpue), mean=logI, sd=s, log=T))
	return(out)
}

line = function(q, alpha, s, sigma, lambda, kappa, omega, w, b){
	k = 1;
	#R0 = 1;
        W0 = w#1/10000;
        B0 = b#10000;
        C0 = 0;
        #
        B = matrix(NA, Y, 1)
        W = matrix(NA, Y, 1)
        R = matrix(NA, Y, 1)
        logR = matrix(NA, Y, 1)
        logI = matrix(NA, Y, 1)
        #
	#//logR[1] = log(alpha) + log(B0-C0) + log(1-beta*gamma*(B0-C0))/gamma; //s
        R[1] = alpha*(B0-C0) #//exp(logR[1]);
        B[1] = R[1] + sigma * (lambda+kappa*W0) * (B0-C0) / W0; 
        W[1] = B[1]/( R[1]/omega + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1] - C[1]/2.);
        for(t in 2:Y){
		#logR[t] = log(alpha) + log(B[t-k]-C[t-k]) + log(1-beta*gamma*(B[t-k]-C[t-k]))/gamma;
                R[t] = alpha*(B[t-k]-C[t-k]) #* (1-beta*(B[t-k]-C[t-k]))#* (1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma) #exp(logR[t]);
                B[t] = R[t] + sigma * (lambda+kappa*W[t-1]) * (B[t-1]-C[t-1]) / W[t-1];
                W[t] = B[t]/( R[t]/omega + sigma*(B[t-1]-C[t-1])/W[t-1] );
                logI[t] = log(q) + log(B[t] - C[t]/2.);
        }
	#
	return( exp(logI) )
}

monitor = function(obj){
	#	
	fitness = na.exclude(obj@fitness)
	where = which(obj@fitness==max(fitness))
	writeLines(sprintf('It=%7d | Best=%3.6e | (%2.3e, %2.3e, %2.3e, %2.3e, %2.3e, %2.3e)', obj@iter, max(fitness), obj@population[where,1], obj@population[where,2], obj@population[where,3], obj@population[where,4], obj@population[where,5], obj@population[where,6]))
	#
	plot(year, cpue)
	lines(year, line(obj@population[where,1], obj@population[where,2], obj@population[where,3], obj@population[where,4], 0.000731, 1.018862, 0.0007, obj@population[where,5], obj@population[where,6]))
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

#loglike     q		  alpha	      beta       s	   sigma       lambda   kappa    omega       W0        B0
#            0.001072705  1.353514    0.4549302  0.9813469 0.1673594
#-4395.19  : 0.0006918396 4.380381    0.3725511  0.9751676 1.06582
#-616.955  : 0.003847164  2.238291    0.4360691  0.9923106 0.1412666 
#-12.70684 : 0.0006445944 0.8109448   0.40096593 0.4205175 0.2504228
#-2.860237 : 0.000489682  0.01069735  0.5843549  0.274137  0.3207339  223.3995
#0.5969832 : 0.0003819466 0.009010002 0.98400676 0.2357711 0.3832153  173.4229 1
#-0.2826202: 0.0002095103 0.02129398  0.5927671  0.2449384 0.06733935 321.2254 8.694229 
#-2.578255 : 0.0001088442 0.09626838  0.6529662  0.2714189 0.05063893 968.5945 17.89035 446.9925
#-1.439356 : 0.0002037883 0.3802532   0.5740516  0.257574  0.01853509 5787.913 32.25317 541.282
#-1.454938 : 0.0002050657 0.3844074 		 0.2594782 0.02478384 4995.887 23.98279 637.7423
#-0.8793316: 0.0001885146 0.1545112   		 0.2514008 0.0540272  7403.512 11.26196 213.8503
#
#1.170625  : 0.0002675379 0.2912155              0.2298821 0.5175151 0.000731  1.018862 0.0007
#2.925287  : 0.000216797  0.1799131              0.2130759 0.6220578 0.000731  1.018862 0.0007
#
#2.873975  : 0.0002250495 0.2895951              0.2135474 0.5093715 0.000731  1.018862 0.0007
#3.149446  : 2.117e-04    2.376e-01              2.112e-01 5.568e-01 0.000731  1.018862 0.0007
#
#4.055288  : 0.0003126413 0.06680737             0.2041795 0.8401128 0.000731  1.018862 0.0007       9974.717  3537.921
#4.056265  : 3.126e-04    6.681e-02              2.028e-01 8.401e-01 0.000731  1.018862 0.0007	     9.881e+04 3.538e+03
#4.056609  : 3.126e-04    6.681e-02              2.035e-01 8.401e-01 0.000731  1.018862 0.0007       9.990e+05 3.538e+03
#4.056731  : 3.125e-04    6.681e-02              2.029e-01 8.401e-01 0.000731  1.018862 0.0007	     4.796e+07 3.538e+03
#2.873975  : 3.116e-04    4.516e-02              1.965e-01 8.919e-01                                 9.998e+07 3.440e+03
sug = rbind(
	c(0.0003126413, 0.06680737, 0.2041795, 0.8401128, 9974.717, 3537.921),	
	c(3.126e-04, 6.681e-02, 2.028e-01, 8.401e-01, 9.881e+04, 3.538e+03),
	c(3.126e-04, 6.681e-02, 2.035e-01, 8.401e-01, 9.990e+05, 3.538e+03)
)
nome = c('q', 'alpha', 's', 'sigma', 'W0', 'B0')
parLow = rep(.Machine$double.eps, dim(sug)[2])
parLow[5] = 1e4
parHi  = c(1e-3, 1e0, 1e0, 1e0, 1e8, 1e4) 
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



