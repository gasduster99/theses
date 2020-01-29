rm(list=ls())

#
library(rootSolve)
library(rstan)
#
source('funcs.r')
#
set.seed(1)

#
#SOLVE
#

#Bs/B0=0.5
#F*/M=1.5 
M  = 0.3#0.15#
Fs = M*1.5
hs = 1-exp(-Fs)
start = c(302, 0.12)
names(start) = c('Cs', 'gamma')
par = data.frame(
	kappa=0,
	wk=0.001,
	wI=0.001,
	Bs=1000,#1500,#
	B0=2000,#3000,#
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
q     = 1
wk    = 0.001
kappa = 0
gaMod = -1
delta = 1-exp(-M)
#
ab    = abGet(hs, Cs, delta, gamma, wk, kappa)
alpha = ab$a
beta  = ab$b
#
Y = 30
C = rep(200, Y)*1000/1500
mod = model(q, alpha, beta, delta, kappa, gamma)
sd = 0.001
cpue = rlnorm(Y, mod$logI, sd)
#
D = list(N=Y, cpue=cpue, C=C, gamma=gaMod, wk=wk, kappa=kappa)

#
#MODEL
#

#
#interface w/ stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
out = stan( file = "numSolve.stan",
        data  = D,
        iter  = 10^4,
        cores = 4,
        init  = function(){list(
                q     = q,
                hStar = hs,
                CStar = Cs,
                delta = delta,
                s     = sd
        )}
)
#
samples = extract(out, permuted=T)
MM = length(samples[[1]])
P = length(samples)-1

#
m = stan_model( file="numSolve.stan" )
optOut = optimizing( m,
        data=D,
        init = function(){list(
                q     = q,
                hStar = hs,
                CStar = Cs,
                delta = delta,
                s     = sd
        )}
)

#
#LOOK
#

#
postR = matrix(NA, nrow=MM, ncol=Y)
postB = matrix(NA, nrow=MM, ncol=Y)
postW = matrix(NA, nrow=MM, ncol=Y)
post  = matrix(NA, nrow=MM, ncol=Y)
pred  = matrix(NA, nrow=MM, ncol=Y)
for(m in 1:MM){
        modOut = model(samples$q[m], samples$alpha[m], samples$beta[m], samples$delta[m], 0, -1)
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

##
#pairs(out)
#

#
pdf(sprintf('./pics/dataB0%d_M%d.pdf', par$B0, M*100))
year=1:Y
plot(year, cpue, ylim=c(min(pred), max(pred)))
polygon(c(year, rev(year)), c(predUpper, rev(predLower)), col='grey', border=NA)
polygon(c(year, rev(year)), c(postUpper, rev(postLower)), col='dimgrey', border=NA)
lines(year, colMeans(pred), lwd=5)
points(year, cpue)
dev.off()

#
#LOOK
#

#
h0 = 0
#
s0 = sGet(h0, samples$delta)
W0 = WGet(s0, wk, kappa)
r0 = rGet(h0, s0, W0, wk)
R0 = RGet(samples$alpha, samples$beta, gaMod, r0)
B0 = BGet(s0, W0, R0, wk)
#
sh = sGet(samples$hStar, samples$delta)
Wh = WGet(sh, wk, kappa)
rh = rGet(samples$hStar, sh, Wh, wk)
Rh = RGet(samples$alpha, samples$beta, gaMod, rh)
Bh = BGet(sh, Wh, Rh, wk)
#
s0Head = sGet(h0, optOut$par['delta'])
W0Head = WGet(s0Head, wk, kappa)
r0Head = rGet(h0, s0Head, W0Head, wk)
R0Head = RGet(optOut$par['alpha'], optOut$par['beta'], gaMod, r0Head)
B0Head = BGet(s0Head, W0Head, R0Head, wk)
#
shHead = sGet(optOut$par['hStar'], optOut$par['delta'])
WhHead = WGet(shHead, wk, kappa)
rhHead = rGet(optOut$par['hStar'], shHead, WhHead, wk)
RhHead = RGet(optOut$par['alpha'], optOut$par['beta'], gaMod, rhHead)
BhHead = BGet(shHead, WhHead, RhHead, wk)
#
FHead = -log(1-optOut$par['hStar'])
MHead = -log(1-optOut$par['delta'])
#
Fmsy = -log(1-samples$hStar)
EM   = -log(1-samples$delta)
#
headx = FHead/MHead   #median(Fmsy/EM)
heady = BhHead/B0Head #median(Bh/B0)
tailx = Fs/M
taily = par$Bs/par$B0
#
pdf(sprintf('./pics/ratioB0%d_M%d.pdf', par$B0, M*100))
plot(Fmsy/EM, Bh/B0, xlim=c(min(Fmsy/EM), max(c(Fs/M, Fmsy/EM))))
arrows(tailx, taily, x1=headx, y1=heady, col='red')
points(tailx, taily, pch=19, col='red')
dev.off()
