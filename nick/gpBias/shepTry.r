rm(list=ls())

#
library(GA)
library(rootSolve)

#
#FUNCTIONS
#

#
dPdt = function(t, P, alpha, beta, gamma, M, catch){
	#
	C = catch[t]
	#
	R = alpha*P/(1+beta*P^(1/gamma))
	out = R - M*P - C #(MM+FF)*P
	#
	return( list(out) )
}

#
shepSRR = function(P, alpha, beta, gamma){ alpha*P/(1+beta*P^(1/gamma)) }

#
PBar = function(ff, alpha, beta, gamma, M){ exp(gamma*log(alpha/(M+ff)-1) - gamma*log(beta)) }
#(alpha/(M+ff)-1)^gamma * beta^-gamma }

#
FMsy = function(alpha, beta, gamma, M){
	#
	FUpper = alpha-M
	out = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))	#
	#uniroot: return(out$root)
	return(out[1])
}

#
getBeta = function(alpha, gamma, M, P0, start=0.01){ 
	optim(start, function(x){ 
			(PBar(0, alpha, x, gamma, M)-P0)^2 
		}, 
		method='Brent', 
		lower=0, 
		upper=10^6
	)$par
}

#
#DATA
#

#
M = 0.2
#delta = 1-exp(-M)
#
P0 = 3000

#               
xi = 2.4	#2          #3.4    #1.1
zeta = 0.35	#0.25     #0.74   #0.6

#
#MODELS
#

#
alpha = 2
gamma = 1
beta = getBeta(alpha, gamma, M, P0)
#
FStar = FMsy(alpha, beta, gamma, M)
PStar = PBar(FStar, alpha, beta, gamma, M)
PZero = PBar(0, alpha, beta, gamma, M)

#
f = function(par, extra){
       #unpack
       alpha = par[1] 
       gamma = par[2]
       #extra
       M = extra$M
       xi = extra$xi
       zeta = extra$zeta
       
       #solve beta against P0
       beta = getBeta(alpha, gamma, M, extra$P0) 

       #par values
       FStar = FMsy(alpha, beta, gamma, M) 
       PStar = PBar(FStar, alpha, beta, gamma, M)
       #ref values
       FSRef = xi*M
       PSRef = zeta*extra$P0

       #
       out = c(FStar-FSRef, PStar-PSRef)
       return( out )
}

#
par = c(alpha, gamma)
extra = data.frame(xi=xi, zeta=zeta, M=M, P0=P0)

#
gaRoot = ga(
	type    = 'real-valued',
        fitness = function(x, extra){ -sum(f(x, extra)^2) }, 
	extra   = extra,
        lower   = c(0, 0),
        upper   = c(5, 5),
        popSize = 100,
        maxiter = 1e6,
        run     = 500,
        optim   = T,
        monitor = T
)
par = gaRoot@solution

#
roots = multiroot(f, par, parms=extra, positive=T)#, maxiter=1e6)

#
alphaR = roots$root[1]
gammaR = roots$root[2]
betaR  = getBeta(alphaR, gammaR, M, P0)
#
FStarR = FMsy(alphaR, betaR, gammaR, M)
PStarR = PBar(FStarR, alphaR, betaR, gammaR, M)
PZeroR = PBar(0, alphaR, betaR, gammaR, M)

#
thresh = 0.05
isNorm = norm(matrix(c(FStarR/M-xi, PStarR/PZeroR-zeta), ncol=2))/norm(matrix(c(xi, zeta), ncol=2))<thresh

##
#fga = function(par, extra){
#	#unpack
#	alpha = par[1]#['alpha']
#	#beta  = #0.003 #par[2]#['beta'] 
#	gamma = par[2]#['gamma']
#	#extra
#	M = extra$M
#	xi = extra$xi
#	zeta = extra$zeta
#	
#	#
#	beta = optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=10^6)$par
#	#beta = opt$par
#
#	#par values
#	#alpha>M+ff
#	FUpper = alpha-M
#	#FStar = uniroot(function(ff){ (ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }, c(0, FUpper))
#	FStar = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper)) 
#	PStar = PBar(FStar, alpha, beta, gamma, M)
#	PZero = PBar(0, alpha, beta, gamma, M)
#	#ref values
#	FSRef = xi*M
#	P0Ref = extra$P0
#	PSRef = zeta*P0Ref	
#
#	#
#	out = c(FStar-FSRef, PStar-PSRef) #, PZero-P0Ref)
#	if( length(out)<2 ){ return(NA) }
#	out = -(sum(out^2))
#	return( out )
#}
#
##
#par = c(alpha, gamma) #data.frame(alpha=alpha, beta=beta, gamma=gamma)
##names(par) = c("alpha", "beta", "gamma")
#extra = data.frame(xi=xi, zeta=zeta, M=M, P0=P0)
#
####
###opt = optim(par, function(p){-f(p, extra)}, 
###	lower = c(M, 0, -10),
###	upper = c(10, 10^4, 10)
###)
##
##
#gaRoot = ga(
#	type    = 'real-valued',
#        fitness = fga, 
#	extra   = extra,
#        lower   = c(0, 0),
#        upper   = c(2, 2),
#        popSize = 100,
#        maxiter = 1e6,
#        run     = 500,
#        optim   = T,
#        monitor = T
#)
#
##
#par = gaRoot@solution
#alpha = par[1]
#gamma = par[2]
#beta = optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=10^6)$par
##
##alpha>M+ff
#FUpper = alpha-M
###(ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }
#FStar = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
#PStar = PBar(FStar, alpha, beta, gamma, M)
#PZero = PBar(0, alpha, beta, gamma, M)
##
#
##
#f = function(par, extra){
#	#unpack
#	alpha = par[1]#['alpha']
#	#beta  = #0.003 #par[2]#['beta'] 
#	gamma = par[2]#['gamma']
#	#extra
#	M = extra$M
#	xi = extra$xi
#	zeta = extra$zeta
#	
#	#
#	beta = optim(beta, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=1000)$par
#	#beta = opt$par
#
#	#par values
#	#alpha>M+ff
#	FUpper = alpha-M
#	#FStar = uniroot(function(ff){ (ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }, c(0, FUpper))
#	FStar = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper)) 
#	PStar = PBar(FStar, alpha, beta, gamma, M)
#	PZero = PBar(0, alpha, beta, gamma, M)
#	#ref values
#	FSRef = xi*M
#	P0Ref = extra$P0
#	PSRef = zeta*P0Ref	
#
#	#
#	out = c(FStar-FSRef, PStar-PSRef) #, PZero-P0Ref)
#	#if( length(out)<2 ){ return(NA) }
#	#out = -(sum(out^2)+1)
#	return( out )
#}
##
#roots = multiroot(f, par, parms=extra, maxiter=1e6, positive=T)
##
#alphaR = roots$root[1]
#gammaR = roots$root[2]
#betaR = optim(0.01, function(x){ (PBar(0, alphaR, x, gammaR, M)-extra$P0)^2 }, method='Brent', lower=0, upper=10^6)$par
##
#Fupper = alphaR-M
#FStarR = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alphaR)+log(gammaR) - (2*log(M+ff)+log(alphaR/(M+ff)-1))) }, c(0, FUpper))
#PStarR = PBar(FStarR, alphaR, betaR, gammaR, M)
#PZeroR = PBar(0, alphaR, betaR, gammaR, M)
##









#
#JUNK
#

##
#abgRootFind = function(){
#	#alpha>M+ff
#	FUpper = alpha-M
#	##(ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }
#	FStar = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper)) 
#	PStar = PBar(FStar, alpha, beta, gamma, M)
#	P0 = PBar(0, alpha, beta, gamma, M)	
#}

##alpha>M+ff
#FUpper = alpha-M
##FStar = uniroot.all(function(ff){ (ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }, c(0, FUpper))
#FStar = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper)) 
#logFStar = uniroot.all(function(lff){ 1 - exp(lff+log(alpha)+log(gamma) - (2*log(M+exp(lff))+log(alpha/(M+exp(lff))-1))) }, c(0, log(FUpper)))
#FStar = exp(logFStar)








