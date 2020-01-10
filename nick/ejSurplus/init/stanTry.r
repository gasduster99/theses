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
	ctch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
)

#
#FIT
#

#interface w/ stan
rstan_options(auto_write = TRUE)
out = stan( file = "stanTry.stan",
	data = D,
	iter = 10^4,
	cores = 4,
	init =  function(){list(
		#logK=8,
		K=2900
		#q=0.00045,
		#r=0.39
	)}
)
#extract samples from stan output
samples = extract(out, permuted=T)
M = length(samples[[1]])
#
B = matrix(NaN, nrow=M, ncol=Y)
meanMean = matrix(NaN, nrow=M, ncol=Y) 
topMean = matrix(NaN, nrow=1, ncol=Y)
botMean = matrix(NaN, nrow=1, ncol=Y)
meanPred = matrix(NaN, nrow=M, ncol=Y)
topPred = matrix(NaN, nrow=1, ncol=Y)
botPred = matrix(NaN, nrow=1, ncol=Y)
#
B[,1] = samples$K 
meanMean[,1] = samples$q*B[,1]
topMean[1] = quantile(meanMean[,1], 0.025)
botMean[1] = quantile(meanMean[,1], 0.975)
meanPred[,1] = rnorm(M, meanMean[,1], samples$s)
topPred[1] = quantile(meanPred[,1], 0.025)
botPred[1] = quantile(meanPred[,1], 0.975)
for(t in 2:Y){
	B[,t] = B[,t-1] + samples$r*B[,t-1]*(1-(B[,t-1]/samples$K)) - D$ctch[t-1]
        meanMean[,t] = samples$q*B[,t]
	topMean[t] = quantile(meanMean[,t], 0.025)
	botMean[t] = quantile(meanMean[,t], 0.975)
	meanPred[,t] = rnorm(M, meanMean[,t], samples$s)
	topPred[t] = quantile(meanPred[,t], 0.025)
        botPred[t] = quantile(meanPred[,t], 0.975)
}

#see fit
plot(year, D$cpue)
polygon(c(year, rev(year)), c(topPred, rev(botPred)), col='grey', border=NA)
polygon(c(year, rev(year)), c(topMean, rev(botMean)), col='dimgrey', border=NA)
lines(year, colMeans(meanPred), lwd=5)
points(year, D$cpue)


