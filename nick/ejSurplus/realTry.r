rm(list=ls())

#
library(rstan)

#
#DATA
#

#
year = 1965:1987; Y = length(year)
#specify data to be imported into stan 
D = list(#use the same R variable names here that will be used in the stan data block
	N = Y,
	cpue = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63),
	#'catch' is an internal variable to stan; no variable names are allowed to be call 'catch'
	C = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
)

#
#FIT
#

#interface w/ stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
out = stan( file = "realTry.stan",
	data = D,
	iter = 10^4,
	cores = 4,
	init = function(){list(
		#K=2900
		q=3.125e-04, #0.000216797, #0.0002675379,
		#sigma=0.5,
		#lambda=3000,
		#kappa=3000,
		alpha=6.681e-02, #0.1799131, #0.2912155,
		#beta=1,
		#gamma=1,
		s=2.029e-01 #0.2130759 #0.2298821
	)}
)

#extract samples from stan output
samples = extract(out, permuted=T)
M = length(samples[[1]])
#
B = matrix(NA, M, Y)
W = matrix(NA, M, Y)
R = matrix(NA, M, Y)
I = matrix(NA, M, Y)
#meanMean = matrix(NaN, nrow=M, ncol=Y) 
topMean = matrix(NaN, nrow=1, ncol=Y)
botMean = matrix(NaN, nrow=1, ncol=Y)
pred = matrix(NaN, nrow=M, ncol=Y)
topPred = matrix(NaN, nrow=1, ncol=Y)
botPred = matrix(NaN, nrow=1, ncol=Y)
#
gamma  = -1;
sigma  = 8.401e-01 #0.5175151
lambda = 0.000731 #1;
kappa  = 1.018862 #1;
omega  = 0.0007;
k = 1;
R0 = 1;
W0 = 4.796e+07 #10000;
B0 = 3.538e+03 #10000;
C0 = 0;
#posterior
R[,1] = samples$alpha*(B0-C0)
B[,1] = R[,1] + sigma * (lambda+kappa*W0) * (B0-C0) / W0;
W[,1] = B[,1]/( R[,1]/omega + sigma*(B0-C0)/W0 );
I[,1] = samples$q * (B[,1] - D$C[1]/2.);
botMean[1] = quantile(I[,1], 0.025)
topMean[1] = quantile(I[,1], 0.975)
#predictive
pred[,1] = rlnorm(M, I[,1], samples$s)
topPred[1] = quantile(pred[,1], 0.975)
botPred[1] = quantile(pred[,1], 0.025)
for(t in 2:Y){
	#posterior
	R[,t] = samples$alpha*(B[,t-k]-D$C[t-k])
        B[,t] = R[,t] + sigma * (lambda+kappa*W[,t-1]) * (B[,t-1]-D$C[t-1]) / W[,t-1]; 
        W[,t] = B[,t]/( R[,t]/omega + sigma*(B[,t-1]-D$C[t-1])/W[t-1] )
        I[,t] = samples$q*(B[,t] - D$C[t]/2.)
	botMean[t] = quantile(I[,t], 0.025)
	topMean[t] = quantile(I[,t], 0.975)
	#predictive
	pred[,t] = rlnorm(M, I[,t], samples$s)
	topPred[t] = quantile(pred[,t], 0.975)
	botPred[t] = quantile(pred[,t], 0.025)
}
#see fit
dev.new()
plot(year, D$cpue, ylim=c(-1, 2))
polygon(c(year, rev(year)), c(topPred, rev(botPred)), col='grey', border=NA)
polygon(c(year, rev(year)), c(topMean, rev(botMean)), col='dimgrey', border=NA)
lines(year, colMeans(I[,1:11]), lwd=5)
points(year, D$cpue)



##
#B[,1] = samples$K 
#meanMean[,1] = samples$q*B[,1]
#topMean[1] = quantile(meanMean[,1], 0.025)
#botMean[1] = quantile(meanMean[,1], 0.975)
#meanPred[,1] = rnorm(M, meanMean[,1], samples$s)
#topPred[1] = quantile(meanPred[,1], 0.025)
#botPred[1] = quantile(meanPred[,1], 0.975)
#for(t in 2:Y){
#	B[,t] = B[,t-1] + samples$r*B[,t-1]*(1-(B[,t-1]/samples$K)) - D$ctch[t-1]
#        meanMean[,t] = samples$q*B[,t]
#	topMean[t] = quantile(meanMean[,t], 0.025)
#	botMean[t] = quantile(meanMean[,t], 0.975)
#	meanPred[,t] = rnorm(M, meanMean[,t], samples$s)
#	topPred[t] = quantile(meanPred[,t], 0.025)
#        botPred[t] = quantile(meanPred[,t], 0.975)
#}



#        #
#        B = matrix(NA, Y, 1)
#        W = matrix(NA, Y, 1)
#        R = matrix(NA, Y, 1)
#        logR = matrix(NA, Y, 1)
#        logI = matrix(NA, Y, 1)
#        #
#        B[1] = R0 + sigma * (lambda+kappa*W0) * (B0-C0) / W0;
#        #//logR[1] = log(alpha) + log(B0-C0) + log(1-beta*gamma*(B0-C0))/gamma; //
#        R[1] = R0; #//exp(logR[1]);
#        W[1] = B[1]/( R[1]/omega + sigma*(B0-C0)/W0 );
#        logI[1] = log(q) + log(B[1] - C[1]/2.);
#        for(t in 2:Y){
#                B[t] = R[t-1] + sigma * (lambda+kappa*W[t-1]) * (B[t-1]-C[t-1]) /
#                #logR[t] = log(alpha) + log(B[t-k]-C[t-k]) + log(1-beta*gamma*(B[t
#                R[t] = alpha*(B[t-k]-C[t-k]) #* (1-beta*(B[t-k]-C[t-k]))#* (1-beta
#                W[t] = B[t]/( R[t]/omega + sigma*(B[t-1]-C[t-1])/W[t-1] );
#                logI[t] = log(q) + log(B[t] - C[t]/2.);
#        }


#see fit
plot(year, D$cpue)
polygon(c(year, rev(year)), c(topPred, rev(botPred)), col='grey', border=NA)
polygon(c(year, rev(year)), c(topMean, rev(botMean)), col='dimgrey', border=NA)
lines(year, colMeans(meanPred), lwd=5)
points(year, D$cpue)


