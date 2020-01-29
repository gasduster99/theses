rm(list=ls())

#
library(rootSolve)
library(rstan)

#
#SOLVE
#

#
f = function(x, pars){
	#
	Cs    = x[1]
	gamma = x[2]
	#
	kappa = pars$kappa
	wk    = pars$wk
	wI    = pars$wI
	Bs    = pars$Bs
	B0    = pars$B0
	M     = pars$M
	hs    = pars$hs
	#
	sigma = exp(-M)
	delta = 1-sigma
	#
	ss    = sigma*(1-hs)
	alpha = (1-ss)/(1-hs) * (1+((gamma*hs)/(1-ss)))^(1/gamma)
	beta  = (hs^2)/((1-hs)*(1-ss+gamma*hs)*Cs)
	#
	h  = 0
	sh = sigma*(1-h)
	Wh = ((1-sh)*wk + (1-kappa)*sh*wI) / (1-kappa*sh)
	rh = (wk/Wh) * ((1-sh)/(1-h))
	Rh = rh/(beta*gamma) * (1-(rh/alpha)^gamma)
	#
	BsEq = (sqrt(alpha)-sqrt(delta)) * (1+sqrt(alpha*delta)-delta) / beta*sqrt(delta) - Bs
	B0Eq = (Wh/wk) * (Rh/(1-sh)) - B0
	#
	#print(c(BsEq, B0Eq))
	return( c(BsEq, B0Eq) )
}

#Bs/B0=0.5
#F*/M=5 
M  = 0.15
Fs = M*5
hs = 1-exp(-Fs)
start = c(270, -1)
names(start) = c('Cs', 'gamma')
par = data.frame(
	kappa=0,
	wk=0.001,
	wI=0.001,
	Bs=1500,
	B0=3000,
	M=M,
	hs=hs
)
#
CG = multiroot(f, start, parms=par,
	maxiter=10^5,
	rtol=.Machine$double.eps,
	atol=.Machine$double.eps,
	ctol=.Machine$double.eps
)$root
Cs = CG[1]
gamma = CG[2]

#
#DATA
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
q = 1
kappa = 0
sigma = exp(-M)
#
ss    = sigma*(1-hs)
alpha = (1-ss)/(1-hs) * (1+((gamma*hs)/(1-ss)))^(1/gamma)
beta  = (hs^2)/((1-hs)*(1-ss+gamma*hs)*Cs)
#
Y=30
C = rep(100, Y)+seq(1,(Y*10), 10)
mod = model(q, alpha, beta, sigma, kappa, gamma)
cpue = rlnorm(Y, mod$logI, 0.001)
#
D = list(N=Y, cpue=cpue, C=C)

#
#MODEL
#

#
#interface w/ stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
out = stan( file = "realTry3.stan",
        data = D,
        iter = 10^4,
        cores = 4,
        init = function(){list(
                q     = q,
                hStar = hs,
                CStar = Cs,
                sigma = sigma,
                s     = 0.001
        )}
)
#
samples = extract(out, permuted=T)
M = length(samples[[1]])
P = length(samples)-1
#
m = stan_model( file="realTry3.stan" )
f = optimizing( m, 
	data=D,  
	init = function(){list(
                q     = q,
                hStar = hs,
                CStar = Cs,
                sigma = sigma,
                s     = 0.001
        )}
)

#
#LOOK
#

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
pairs(out)
#
dev.new()
year=1:Y
plot(year, cpue, ylim=c(min(pred), max(pred)))
polygon(c(year, rev(year)), c(predUpper, rev(predLower)), col='grey', border=NA)
polygon(c(year, rev(year)), c(postUpper, rev(postLower)), col='dimgrey', border=NA)
lines(year, colMeans(pred), lwd=5)
points(year, cpue)

#
wk = 0.001
gamma = -1
kappa = 0
#
h  = 0
sh = samples$sigma*(1-h)
W0 = ((1-sh)*wk + (1-kappa)*sh*wk)/(1-kappa*samples$sigma)
rh = (wk/W0)*(1-sh)/(1-h)
R0 = (1-(rh/samples$alpha)^gamma)*rh/(samples$beta*gamma) #(1-exp(gamma*(log(rh)-log(alpha))))*(rh/(beta*
B0 = (W0/wk)*(R0/(1-sh))
#
rStar = sqrt(samples$alpha*(1-samples$sigma))
hStar = rStar - (1-samples$sigma)
hStar = hStar/(1+hStar)
CStar = ((sqrt(samples$alpha) - sqrt(1-samples$sigma))^2) / (samples$beta)
RStar = ((sqrt(samples$alpha) - sqrt(1-samples$sigma))    / (samples$beta)) * sqrt(samples$alpha)
SStar =  (sqrt(samples$alpha) - sqrt(1-samples$sigma))    / (samples$beta*sqrt(1-samples$sigma))
BStar = SStar * (1+sqrt(samples$alpha*(1-samples$sigma))-(1-samples$sigma))
#
Fmsy = -log(1-hStar)
EM = -log(samples$sigma)
#
dev.new()
plot(Fmsy/EM, BStar/B0)

#
wk = 0.001
gamma = -1
kappa = 0
#
hStar = summary(out)$summary[3,6] #f$par['hStar']
sigma = summary(out)$summary[4,6] #f$par['sigma']
alpha = summary(out)$summary[7,6] #f$par['alpha']
beta  = summary(out)$summary[8,6] #f$par['beta']
#
h  = 0
sh = sigma*(1-h)
W0 = ((1-sh)*wk + (1-kappa)*sh*wk)/(1-kappa*sigma)
rh = (wk/W0)*(1-sh)/(1-h)
R0 = (1-(rh/alpha)^gamma)*rh/(beta*gamma) #(1-exp(gamma*(log(rh)-log(alpha))))*(rh/(beta*
B0 = (W0/wk)*(R0/(1-sh))
#
rStar = sqrt(alpha*(1-sigma))
hStar = rStar - (1-sigma)
hStar = hStar/(1+hStar)
CStar = ((sqrt(alpha) - sqrt(1-sigma))^2) / (beta)
RStar = ((sqrt(alpha) - sqrt(1-sigma))    / (beta)) * sqrt(alpha)
SStar =  (sqrt(alpha) - sqrt(1-sigma))    / (beta*sqrt(1-sigma))
BStar = SStar * (1+sqrt(alpha*(1-sigma))-(1-sigma))
#
Fmsy = -log(1-hStar)
EM = -log(sigma)
points(Fs/0.15, par$Bs/par$B0, col='red')
arrows(Fs/0.15, par$Bs/par$B0, x1=Fmsy/EM, y1=BStar/B0, col='red')


