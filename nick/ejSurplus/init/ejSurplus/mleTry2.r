rm(list=ls())

#
#FUNCTIONS
#

BNext = function(BNow, r, K, CNow){
	#BNow + r*BNow*(1-log(BNow/K)) - CNow
	one = BNow 
	two1 = r*BNow
	two2 = 1-exp(log(BNow)-log(K))
	thr = -CNow
	#
	out = one + two1*two2 + thr
	return(out)
}

#
ll = function(par, x){
	#parameters
	q = par[1]
	r = par[2]
	K = par[3]
	s = par[4]
	B0= par[5]
	#
	#B0 = K
	C0 = 0#C[1]
	#
	Bs = matrix(NaN, nrow=Y, ncol=1)
	Bs[1] = BNext(B0, r, K, C0)
	for(t in 1:Y){ Bs[t+1]=BNext(Bs[t], r, K, C[t]) }
	#Is = q*Bs
	logIs = log(q) + log(Bs)
	#out = -sum(dnorm(log(y), mean=log(Is), sd=s, log=T))
	out = -sum(dnorm(log(x), mean=logIs, sd=s, log=T))
	#print(c(q, r, K, s, B0))
	#print(logIs)
	#print('')
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
#optimization
s0 = 0.1
q0 = 0.001
r0 = 0.1
K0 = 10000
B0 = 10000
par    = c(q0, r0, K0, s0, B0)
parLow = rep(.Machine$double.eps, 5)
parHi  = c(1e0, 1e0, 1e6, 1e0, 1e6)
optimOut = optim(par, ll, x=y,
	method  = "L-BFGS-B",
	lower   = parLow,
	upper   = parHi,
	control = list(
		maxit=1e9,
		ndeps=rep(.Machine$double.eps, 5),
		pgtol=.Machine$double.eps,
		fnscale=.Machine$double.eps
	)
)
#output
sHat  = 0.12   	#optimOut$par[4]	#
qHat  = 0.00045	#optimOut$par[1]	#
rHat  = 0.39   	#optimOut$par[2]	#
KHat  = 2709   	#optimOut$par[3]	#
B0Hat = KHat	#optimOut$par[5] 	#
#
writeLines(sprintf('rHat: %f', rHat))
writeLines(sprintf('KHat: %f', KHat))
writeLines(sprintf('qHat: %f', qHat))
writeLines(sprintf('sHat: %f', sHat))
writeLines(sprintf('B0Hat: %f', B0Hat))
#
IHat = matrix(NaN, nrow=1, ncol=Y); colnames(IHat) = year
BHat = matrix(NaN, nrow=1, ncol=Y); colnames(BHat) = year
BHat[1] = BNext(B0Hat, rHat, KHat, 0)
for(t in 1:(Y-1)){
        BHat[t+1] = BNext(BHat[t], rHat, KHat, C[t])
}
IHat = qHat*BHat
#
plot(year, y)
lines(year, IHat)

