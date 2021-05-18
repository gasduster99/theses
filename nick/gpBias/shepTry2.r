rm(list=ls())

#
library(GA)
library(rootSolve)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
strongRoot = function(f, par, extra, howGood, lower=c(0, 0), upper=c(2, 2), monitor=F){
	#f	: a function computing a system of equations to find a root of
	#par	: a vactor of starting values for computing the roots
	#extra	: a data.frame of other incidental values required for f
	#howGood: a function evaluating how good each root found is
	#lower	: bounds of root search
	#upper	: bounds of root search
	#monitor: do you want to watch ga progress 
	#
	#value	: a root of f as a vector

	#
	wrap = function(x, extra){ 
		#send unwanted output to /dev/null
		capture.output(
		#appearantly you need the <- inside capture.output else out will not be defined outside
		out <- tryCatch({
			#parMR = multiroot(f, x, parms=extra, maxiter=1e6, positive=T)$root
			out = howGood(x, extra)
		}, error=function(err){
			out = -Inf
		}), file="/dev/null") 
		#
		return(out)
	}

	#
	dm = length(lower)
	vol = prod(upper-lower)
	gaRoot = ga(
	        type    = 'real-valued',
	        fitness = wrap,
	        extra   = extra,
	        lower   = lower,
	        upper   = upper,
	        popSize = 100*dm,
	        maxiter = 1e6,
	        run     = 200*dm,
		optim   = T,
	        monitor = monitor#,suggestions = t(par)
	)
	parGA = gaRoot@solution[1,]
	isGAGood = isGood(parGA, extra)
	howGAGood = howGood(parGA, extra)
	#print(isGood(parGA, extra))
	#print(howGood(parGA, extra))
	#print(parGA)
	##
	parMR = multiroot(f, parGA, parms=extra, maxiter=1e6, positive=T)$root
	isMRGood = isGood(parMR, extra)
	howMRGood = howGood(parMR, extra)
	#print( isGood(parMR, extra) )	
	#print( howGood(parMR, extra) )
	#print( parMR )
	
	#
	l = c(GA=howGAGood, MR=howMRGood)
	who = which(l==min(l, na.rm=T))
	
	#
	if( isGAGood & isMRGood ){ return(c(GA=parGA, MR=parMR)[who]) }
	if( isGAGood & !isMRGood ){ return(c(GA=parGA)) }
	if( !isGAGood & isMRGood ){ return(c(MR=parMR)) }

	#
	return(c(GA=parGA, MR=parMR)[who])
}


#
dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
        #
        C = catch[t]
        #
        R = exp(lalpha)*P/(1+exp(lbeta)*P^(1/gamma))
        out = R - M*P - C
        #
        return( list(out) )
}
##
#dPdt = function(t, P, alpha, beta, gamma, M, catch){
#	#
#	C = catch[t]
#	#
#	R = alpha*P/(1+beta*P^(1/gamma))
#	out = R - M*P - C 
#	#
#	return( list(out) )
#}

#
shepSRR = function(P, alpha, beta, gamma){ alpha*P/(1+beta*P^(1/gamma)) }

#
#P0Funk = function(alpha, beta, gamma, M){ (alpha/(M)-1)^gamma * beta^-gamma }
#exp(gamma*log(alpha/(M+ff)-1) - gamma*log(beta)) }
PBar = function(ff, alpha, beta, gamma, M){ (alpha/(M+ff)-1)^gamma * beta^-gamma }

#
getBeta = function(alpha, gamma, M, P0){
	optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-P0)^2 }, method='Brent', lower=0, upper=10^6)$par
}

#
FMsy = function(alpha, gamma, M){
	#
	FUpper = alpha-M
	root = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
	#
	if(length(root)<1){ return(NA) }
	if(length(root)>1){ return(root[1]) }
	#
	return(root)
}

#
f = function(par, extra){	
	#unpack
	alpha = par[1]
	gamma = par[2]
	#extra
	M = extra$M
	xi = extra$xi
	zeta = extra$zeta
	
	#	
	beta = getBeta(alpha, gamma, M, extra$P0) 
	
	#par values 
	FStar = FMsy(alpha, gamma, M) 
	PStar = PBar(FStar, alpha, beta, gamma, M)
	#ref values
	FSRef = xi*M
	PSRef = zeta*extra$P0	

	#
	out = c(FStar-FSRef, PStar-PSRef)	 
	return( out )
}

#
howGood = function(par, extra){
	#unpack
        alpha = par[1]
        gamma = par[2]
	beta = getBeta(alpha, gamma, M, extra$P0)
	#compute
	PZero = PBar(0, alpha, beta, gamma, M)
	fOut = f(par, extra)
	#handel case of some numerical issue in either PBar or FMsy
	if( length(PZero)<1 ){ return(-Inf) } 
	if( length(fOut)<2 ){ return(-Inf) }
	#c(FStar/M-xi, PStar/PZero-zeta)
	refComp = c(fOut[1]/extra$M, (fOut[2]+extra$P0*extra$zeta)/PZero-extra$zeta)
	propNorm = norm(matrix(refComp, ncol=2))/norm(matrix(c(extra$xi, extra$zeta), ncol=2))
	return( -propNorm )
}
isGood = function(par, extra, thresh=0.05){ 
	#
	out = (-howGood(par, extra))<thresh 
	if(is.na(out)){ out=F }
	#
	return( out )
}

#
#DATA
#

#
hake  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(hake)

#
M = 0.2
#delta = 1-exp(-M)
#
P0 = 3000

#               
xi = 3.4		#2          #3.4    #1.1
zeta = 0.75	#0.25     #0.74   #0.6

#
#ROOTS
#

#
alpha = 2
gamma = 1
beta  = getBeta(alpha, gamma, M, P0) #0.003
#alpha>M+ff
FUpper = alpha-M
FStar = FMsy(alpha, gamma, M) #uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper)) 
PStar = PBar(FStar, alpha, beta, gamma, M)
PZero = PBar(0, alpha, beta, gamma, M)

#something about GA requires par to index be indexed wiout names par could not be a data.frame
#data.frame(alpha=alpha, beta=beta, gamma=gamma)
#names(par) = c("alpha", "beta", "gamma")
par = c(alpha, gamma) 
extra = data.frame(xi=xi, zeta=zeta, M=M, P0=P0)

#
MR = multiroot(f, par, parms=extra, maxiter=1e6, positive=T)
par = MR$root
#
if( !isGood(par, extra) ){
	par = strongRoot(f, par, extra, howGood, lower=c(0, 0), upper=c(5, 5), monitor=T)
}
#
writeLines('\nRoots')
print(isGood(par, extra))
print(howGood(par, extra))

#
#POP DYN
#

#
alpha = par[1]
gamma = par[2]
beta = getBeta(alpha, gamma, M, P0)

#
truth = prodModel$new(
	dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(0, exp(lalpha), exp(lbeta), gamma, M)}, #PBar(0, alpha, beta, gamma, M)}, #model
        time=1:TT, catch=catch, M=M,			#constants
	lalpha=log(alpha), lbeta=log(beta), gamma=gamma,#parameters
        lq=log(0.004), lsdo=log(0.1160256),			#nuisance parameters
        xi=xi, zeta=zeta				#other incidentals to carry along
)
truth$iterate()

#cpue 
dat = rlnorm(TT, truth$lq+log(truth$N), exp(truth$lsdo))

##
#layout(matrix(1:2, ncol=2))
#curve(shepSRR(x, truth$alpha, truth$beta, truth$gamma), 0, 20000)
#truth$plotMean()

#
fit = prodModel$new(
        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(0, exp(lalpha), exp(lbeta), gamma, M)}, #PBar(0, alpha, beta, gamma, M)}, #model
        time=1:TT, catch=catch, M=M,  			#constants
        lalpha=log(alpha), lbeta=log(beta), gamma=1,	#parameters
        lq=log(0.004), lsdo=log(0.1160256),			#nuisance parameters
        xi=xi, zeta=zeta				#other incidentals to carry along
)
#
optAns = fit$optimize(dat,
        c('lsdo', 'lalpha', 'lbeta'), #lq is optimized via profile likelihood internally
        lower   = c(log(0.001), log(M), log(10^-6)),
        upper   = c(log(1), log(100), log(10)),
        gaBoost = list(run=10, parallel=FALSE, popSize=10^3),
	persistFor = 5
)



#
#JUNK
#

##
#fga = function(par, extra){
#	#unpack #NOTE: names didn't work, but numbered index does work
#	alpha = par[1]#['alpha'] 
#	gamma = par[2]#['gamma']
#	#extra
#	M = extra$M
#	xi = extra$xi
#	zeta = extra$zeta
#	
#	#
#	beta = optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=10^6)$par
#
#	#par values
#	#alpha>M+ff
#	FUpper = alpha-M
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
#FStarG = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
#PStarG = PBar(FStarG, alpha, beta, gamma, M)
#PZeroG = PBar(0, alpha, beta, gamma, M)
##
#thresh = 0.05
#isNormG = norm(matrix(c(FStarG/M-xi, PStarG/PZeroG-zeta), ncol=2))/norm(matrix(c(xi, zeta), ncol=2))<thresh
#
##
#f = function(par, extra){
#	#unpack
#	alpha = par[1]#['alpha']
#	gamma = par[2]#['gamma']
#	#extra
#	M = extra$M
#	xi = extra$xi
#	zeta = extra$zeta
#	
#	#
#	beta = optim(beta, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=1000)$par
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
#	out = c(FStar-FSRef, PStar-PSRef) 
#	return( out )
#}
#
#
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
#thresh = 0.05
#isNormR = norm(matrix(c(FStarR/M-xi, PStarR/PZeroR-zeta), ncol=2))/norm(matrix(c(xi, zeta), ncol=2))<thresh#
##

#print( f(par, extra) )
#print( -sum(f(par, extra)^2) )

##
#gaRoot = ga(
#         type    = 'real-valued',
#         fitness = function(x, extra){ -sum(f(x, extra)^2) },
#         extra   = extra,
#         lower   = c(0, 0),
#         upper   = c(2, 2),
#         popSize = 200,
#         maxiter = 1e6,
#         run     = 500,
#         optim   = T,
#         monitor = T
#)

#
#out = tryCatch({
#}, warning=function(war){
#		out = c(0, 0)
#	}, error=function(err){
#		out = c(0, 0)
#	})


	##
	#alpha = par[1]
	#gamma = par[2]
	##
	#M = extra$M
	#xi = extra$xi
	#zeta = extra$zeta
	#
	##
	##optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=10^6)$par
	#beta = getBeta(alpha, gamma, M, extra$P0) 
	#
	##
	##alpha>M+ff
	##FUpper = alpha-M
	###(ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }
	##uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
	#FStar = FMsy(alpha, gamma, M) 
	#PStar = PBar(FStar, alpha, beta, gamma, M)
	#PZero = PBar(0, alpha, beta, gamma, M)
	##
	#if( length(c(FStar, PStar, PZero))<3 ){ return(FALSE) }


#alpha>M+ff
	#FUpper = alpha-M
	#FStar = uniroot(function(ff){ (ff*alpha*gamma)/((M+ff)^2*(alpha/(M+ff)-1)) }, c(0, FUpper))
	#uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))


#optim(beta, function(x){ (PBar(0, alpha, x, gamma, M)-extra$P0)^2 }, method='Brent', lower=0, upper=1000)$par

##
#wrap = function(x, extra){
#        #
#	out = tryCatch({
#                parMR = multiroot(f, x, parms=extra, maxiter=1e6, positive=T)$root
#		out = howGood(parMR, extra)
#        }, error=function(err){
#                print(err)
#		out = -Inf
#        })
#        #
#        return(out)
#}
