rm(list=ls())

#
library(rstan)

#
#FUNCTIONS
#

#
model = function(q, alpha, beta, sigma, kappa, gamma){
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
        B[1] = R[1] + (((1-kappa)*wI + kappa*W0)/W0) * sigma*(B0-C0) #exp(log( (1-kappa)*wI+kappa*W0 ) + log(sigma) + log(B0) - log(W0
        W[1] = B[1]/( R[1]/wk + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1] - C[1]);
        for(t in 2:Y){
                if(t>k){ R[t] = alpha*(B[t-k]-C[t-k])*(1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma)
                }else{   R[t] = R0 }
                B[t] = R[t] + (((1-kappa)*wI + kappa*W[t-1])/W[t-1]) * sigma*(B[t-1]-C[t-1]) #exp(log( (1-kappa)*
                W[t] = B[t]/( R[t]/wk + sigma*(B[t-1]-C[t-1])/W[t-1] );
                logI[t] = log(q) + log(B[t] - C[t]);
        }
        #
	out = list(
		R=R,
		B=B,
		W=W,
		logI=logI
	)
	#
        return( out )
}

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
out = stan( file = "realTry2.stan",
        data = D,
        iter = 10^4,
        cores = 4,
        init = function(){list(
                q     = 4.666e-04, 
                alpha = 3.012e+00, 
                beta  = 6.530e-03,
		sigma = 8.491e-01,
                s     = 1.175e-01 
        )}
)
#
samples = extract(out, permuted=T)
M = length(samples[[1]])
P = length(samples)-1

#
#LOOK
#

#
C = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
#
postR = matrix(NA, nrow=M, ncol=Y)
postB = matrix(NA, nrow=M, ncol=Y)
postW = matrix(NA, nrow=M, ncol=Y)
post = matrix(NA, nrow=M, ncol=Y)
pred = matrix(NA, nrow=M, ncol=Y)
for(m in 1:M){
	modOut = model(samples$q[m], samples$alpha[m], samples$beta[m], samples$sigma[m], 0, -1)
	postR[m,] = modOut$R
	postB[m,] = modOut$B
	postW[m,] = modOut$W
	post[m,] = exp(modOut$logI)	
	pred[m,] = rlnorm(Y, modOut$logI, samples$s[m])
}
#
postUpper = matrix(NA, nrow=1, ncol=Y)
postLower = matrix(NA, nrow=1, ncol=Y)
#
predUpper = matrix(NA, nrow=1, ncol=Y)
predLower = matrix(NA, nrow=1, ncol=Y)
for(t in 1:Y){
	#
	postUpper[t] = quantile(post[,t], 0.025)
	postLower[t] = quantile(post[,t], 0.975)
	#
	predUpper[t] = quantile(pred[,t], 0.025)
	predLower[t] = quantile(pred[,t], 0.975)
}

#
dev.new()
plot(year, D$cpue, ylim=c(0, 2))
polygon(c(year, rev(year)), c(predUpper, rev(predLower)), col='grey', border=NA)
polygon(c(year, rev(year)), c(postUpper, rev(postLower)), col='dimgrey', border=NA)
lines(year, colMeans(pred), lwd=5)
points(year, D$cpue)

#
save.image('save.RData')

