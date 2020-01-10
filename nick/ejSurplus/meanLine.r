line = function(q, alpha, s, sigma, lambda, kappa, omega){
	k = 1;
	R0 = 1;
        W0 = 10000;
        B0 = 10000;
        C0 = 0;
        #
        B = matrix(NA, Y, 1)
        W = matrix(NA, Y, 1)
        R = matrix(NA, Y, 1)
        logR = matrix(NA, Y, 1)
        logI = matrix(NA, Y, 1)
        #
        B[1] = R0 + sigma * (lambda+kappa*W0) * (B0-C0) / W0;
        #//logR[1] = log(alpha) + log(B0-C0) + log(1-beta*gamma*(B0-C0))/gamma; //s
        R[1] = R0; #//exp(logR[1]);
        W[1] = B[1]/( R[1]/omega + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1] - C[1]/2.);
        for(t in 2:Y){
                B[t] = R[t-1] + sigma * (lambda+kappa*W[t-1]) * (B[t-1]-C[t-1]) / W[t-1];
                #logR[t] = log(alpha) + log(B[t-k]-C[t-k]) + log(1-beta*gamma*(B[t-k]-C[t-k]))/gamma;
                R[t] = alpha*(B[t-k]-C[t-k]) #* (1-beta*(B[t-k]-C[t-k]))#* (1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma) #exp(logR[t]);
                W[t] = B[t]/( R[t]/omega + sigma*(B[t-1]-C[t-1])/W[t-1] );
                logI[t] = log(q) + log(B[t] - C[t]/2.);
        }
	#
	return( exp(logI) )
}

#
year = 1965:1987;
Y = length(year)
cpue = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
C = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
#
plot(year, cpue)
lines(year, line(0.0001885146, 0.1545112, 0.2514008, 0.0540272, 7403.512, 11.26196, 213.8503), col='red')
lines(year, line(0.000216797, 0.1799131, 0.2130759, 0.6220578, 0.000731, 1.018862, 0.0007))
