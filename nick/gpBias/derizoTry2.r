rm(list=ls())

#
library(GA)
library(rootSolve)

#
#FUNCTIONS
#

#
strongRoot = function(f, par, extra, howGood, lower=c(0, 0), upper=c(2, 2), monitor=F){
        #f      : a function computing a system of equations to find a root of
        #par    : a vactor of starting values for computing the roots
        #extra  : a data.frame of other incidental values required for f
        #howGood: a function evaluating how good each root found is
        #lower  : bounds of root search
        #upper  : bounds of root search
        #monitor: do you want to watch ga progress 
        #
        #value  : a root of f as a vector

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
        capture.output( parMR <- suppressWarnings(multiroot(f, parGA, parms=extra, maxiter=1e6, positive=T)$root), file="/dev/null")
        isMRGood = isGood(parMR, extra)
        howMRGood = howGood(parMR, extra)
        #print( isGood(parMR, extra) )  
        #print( howGood(parMR, extra) )
        #print( parMR )

        #
        l = c(GA=howGAGood, MR=howMRGood)
        who = which(l==min(l, na.rm=T))

        #
        if( isGAGood & isMRGood ){ return(list(GA=parGA, MR=parMR)[[who[1]]]) }
        if( isGAGood & !isMRGood ){ return(parGA) }
        if( !isGAGood & isMRGood ){ return(parMR) }

        #	
        return(list(GA=parGA, MR=parMR)[[who[1]]])
}


#
dPdt = function(t, P, alpha, beta, gamma, M, catch){
	#
	C = catch[t]
	#
	R = alpha*P/(1+beta*P)^(gamma)
	out = R - M*P - C #(MM+FF)*P
	#
	return( list(out) )
}

#
derizoSRR = function(P, alpha, beta, gamma){ alpha*P/(1+beta*P)^(gamma) }

#
PBar = function(ff, alpha, beta, gamma, M){ ((alpha/(M+ff))^(1/gamma)-1) * 1/beta }

#
FMsy = function(alpha, beta, gamma, M){
	#
	#FUpper = alpha-M
	out = uniroot.all(function(ff){ 1 - (alpha/(M+ff))^(-1/gamma) - (ff/gamma)*(1/(M+ff)) }, c(0, 10))
	#
	return(out)	
}


#
f = function(par, extra){
       	#unpack
       	alpha = par[1] 
       	beta  = par[2]
	gamma = par[3]
       	#extra
       	M = extra$M
       	xi = extra$xi
       	zeta = extra$zeta
       	
       	#par values
       	FStar = FMsy(alpha, beta, gamma, M) 
       	PStar = PBar(FStar, alpha, beta, gamma, M)
	PZero = PBar(0, alpha, beta, gamma, M)
	#ref values
       	FSRef = xi*M
       	PSRef = zeta*extra$P0
	PZRef = extra$P0

       	#
       	out = c(FStar-FSRef, PStar-PSRef, PZero-PZRef)
       	return( out )
}

#
howGood = function(par, extra){
        #unpack
        alpha = par[1]
	beta  = par[2]
        gamma = par[3]
        #compute
        fOut = f(par, extra)
	#
	PStar = fOut[2]+extra$P0*extra$zeta
	PZero = fOut[3]+extra$P0
        #handel case of some numerical issue in either PBar or FMsy 
        if( length(fOut)<3 ){ return(-Inf) }
	if( PZero<0 ){ return(-Inf) }
	#c(FStar/M-xi, PStar/PZero-zeta)
        refComp = c(fOut[1]/extra$M, PStar/PZero-extra$zeta)
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
M = 0.2
#delta = 1-exp(-M)
#
P0 = 3000

#               
xi = 2		#2          #3.4    #1.1
zeta = 0.25	#0.25     #0.74   #0.6

#
#MODELS
#

#
alpha = 1
gamma = 1
beta = 0.001 #getBeta(alpha, gamma, M, P0)
#
FStar = FMsy(alpha, beta, gamma, M)
PStar = PBar(FStar, alpha, beta, gamma, M)
PZero = PBar(0, alpha, beta, gamma, M)

#
par = c(alpha, beta, gamma)
extra = data.frame(xi=xi, zeta=zeta, M=M, P0=P0)

#
sr = strongRoot(f, par, extra, howGood, lower=c(0, 0, -5), upper=c(10, 10, 10), monitor=T)
alpha = sr[1]
beta  = sr[2]
gamma = sr[3]
#
FStar = FMsy(alpha, beta, gamma, M)
PStar = PBar(FStar, alpha, beta, gamma, M)
PZero = PBar(0, alpha, beta, gamma, M)
#
print( c(xi, zeta) )
print( c(FStar/M, PStar/PZero) )
print( c(alpha, beta, gamma) )
print( howGood(sr, extra) )
print( isGood(sr, extra, thresh=0.01) )


##
#lower   = c(0, 0, 0)
#upper   = c(100, 5, 5)
#vol = prod(upper-lower)
#gaRoot = ga(
#	type    = 'real-valued',
#        fitness = function(x, extra){ -sum(f(x, extra)^2) }, 
#	extra   = extra,
#        lower   = lower,
#        upper   = upper,
#        popSize = 2000, #*vol,
#        maxiter = 1e6,
#        run     = 500,
#        optim   = T,
#        monitor = T
#)
#par = gaRoot@solution
#
##
#alphaG = par[1] #roots$root[1]
#gammaG = par[3] #roots$root[2]
#betaG  = par[2] #getBeta(alphaR, gammaR, M, P0)
##
#FStarG = FMsy(alphaG, betaG, gammaG, M)
#PStarG = PBar(FStarG, alphaG, betaG, gammaG, M)
#PZeroG = PBar(0, alphaG, betaG, gammaG, M)
#
##
#thresh = 0.05
#isNormG = norm(matrix(c(FStarG/M-xi, PStarG/PZeroG-zeta), ncol=2))/norm(matrix(c(xi, zeta), ncol=2))<thresh
#
##
#roots = multiroot(f, par, parms=extra, positive=T)#, maxiter=1e6)
#
##
#alphaR = roots$root[1]
#gammaR = roots$root[3]
#betaR  = roots$root[2]
##
#FStarR = FMsy(alphaR, betaR, gammaR, M)
#PStarR = PBar(FStarR, alphaR, betaR, gammaR, M)
#PZeroR = PBar(0, alphaR, betaR, gammaR, M)
#
##
#thresh = 0.05
#isNormR = norm(matrix(c(FStarR/M-xi, PStarR/PZeroR-zeta), ncol=2))/norm(matrix(c(xi, zeta), ncol=2))<thresh









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

#out = optim(2*M, function(ff){ 
	#		-exp(ff)*PBar(exp(ff), alpha, beta, gamma, M) 
	#	},
	#	#gr = function(ff){
	#	#	-PBar(exp())
	#	#}
	#	method='Brent',
	#	lower=-10,
	#	upper=10 #FUpper 
	#)
	#uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))	#
	#uniroot: return(out$root)
	#return(exp(out$par))


##
#getBeta = function(alpha, gamma, M, P0, start=0.01){ 
#	optim(start, function(x){ 
#			(PBar(0, alpha, x, gamma, M)-P0)^2 
#		}, 
#		method='Brent', 
#		lower=0, 
#		upper=10^6
#	)$par
#}


##
#lower   = c(M, 0, 0)
#upper   = c(10, 5, 5)
#jdeRoot = JDEoptim(lower, upper, 
#	function(x, extra){ -sum(f(x, extra)^2) },
#	extra = extra
#)

##
#f = function(par, extra){
#       	#unpack
#       	alpha = par[1] 
#       	beta  = par[2]
#	gamma = par[3]
#       	#extra
#       	M = extra$M
#       	xi = extra$xi
#       	zeta = extra$zeta
#       	
#       	#par values
#       	FStar = FMsy(alpha, beta, gamma, M) 
#       	PStar = PBar(FStar, alpha, beta, gamma, M)
#	PZero = PBar(0, alpha, beta, gamma, M)
#	#ref values
#       	FSRef = xi*M
#       	PSRef = zeta*extra$P0
#	PZRef = extra$P0
#
#       	#
#       	out = c(FStar-FSRef, PStar-PSRef, PZero-PZRef)
#       	return( out )
#}

#PZero = PBar(0, alpha, beta, gamma, M)
	#if( length(PZero)<1 ){ return(-Inf) }

