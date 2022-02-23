rm(list=ls())

#
library(GA)
library(numDeriv)
library(rootSolve)

#
#FUNCTIONS
#

#
strongRoot = function(f, par, extra, howGood, lower=rep(0, length(par)), upper=rep(2, length(par)), monitor=F, gaOnly=F, gaOptim=T){
        #f      : a function computing a system of equations to find a root of
        #par    : a vactor of starting values for computing the roots
        #extra  : a data.frame of other incidental values required for f
        #howGood: a function evaluating how good each root found is
        #lower  : bounds of root search
        #upper  : bounds of root search
        #monitor: do you want to watch ga progress 
        #gaOnly : a hook to skip the multiroot try
        #gaOptim: a hook to try gradient descent attempts
        #
        #value  : a root of f as a vector

        #
        wrap = function(x, extra){
                #send unwanted output to /dev/null
                capture.output(
                #appearantly you need the <- inside capture.output else out will not be defined outside
                out <- tryCatch({
                        #parMR = multiroot(f, x, parms=extra, maxiter=1e6, positive=T)$root
                        out = -sum(f(x)^2) #howGood(x, extra)
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
                optim   = gaOptim,
                monitor = monitor#,suggestions = t(par)
        )
        parGA = gaRoot@solution[1,]
        #isGAGood = isGood(parGA, extra)
        #howGAGood = howGood(parGA, extra)

        ##
        #parMR = NA
        #howMRGood = NA
        #isMRGood = F
        if( !gaOnly ){
                #
                parMR = multiroot(f, parGA, parms=extra, maxiter=1e6, positive=T)$root
                ##
                #isMRGood = isGood(parMR, extra)
                #howMRGood = howGood(parMR, extra)
        }

	#
        l = c(wrap(parGA), wrap(parMR))
	who = which(l==max(l, na.rm=T))
	
	#
	return(l[who])

	##
        #l = c(GA=howGAGood, MR=howMRGood)
        #who = which(l==max(l, na.rm=T))

        ##
        #if( isGAGood & isMRGood ){ return(list(GA=parGA, MR=parMR)[[who]]) }
        #if( isGAGood & !isMRGood ){ return(c(GA=parGA)) }
        #if( !isGAGood & isMRGood ){ return(c(MR=parMR)) }

        ##both are not good, but return best try anyway
        #return(list(GA=parGA, MR=parMR)[[who]])
}

#
SRR = function(B, alpha, beta, gamma){
	alpha*B/(1+beta*B)^gamma
}

#
PBar = function(alpha, beta, gamma, ff, M){ ((alpha/(M+ff))^(1/gamma)-1)/beta }

#
getBeta = function(alpha, gamma, M, B0){
	(1-(alpha/M)^(1/gamma))/B0/gamma
}
b = Vectorize(getBeta, c("alpha", "gamma"))

#
getAlpha = function(gamma, ff, M){
	(M+ff)*(1-(ff/gamma)*(1/(M+ff)))^(-gamma)	
	#(M+ff)*(1-(ff/(M+ff)/gamma))^(-gamma)
}
a = Vectorize(getAlpha, "gamma")

#
getGamma = function(zeta, ff, M){
	uniroot.all(function(gamma){ ((ff/(M+ff))*(1-zeta)/(zeta*gamma)+1)^gamma - (M+ff)/M }, c(-1, 100))
	#strongRoot()
}

#
getZeta = function(gamma, ff, M){
	(ff/(M+ff))/( gamma*((ff+M)/(M))^(1/gamma) - gamma + (ff/(M+ff)) )
}
##
#getZeta = function(gamma){#, ff, M){
#        (ff/(mm+ff))/( gamma*((ff+mm)/(mm))^(1/gamma) - gamma + (ff/(mm+ff)) )
#}
z = Vectorize(getZeta, "gamma")


#
FMsy = function(alpha, gamma, M){
        #
        FUpper = 1 #alpha-M
        root = uniroot.all(function(ff){ 1 - ff/gamma/(M+ff) - (alpha/(M+ff))^(-1/gamma) }, c(0, FUpper))
        #
        if(length(root)<1){ return(NA) }
        if(length(root)>1){ return(root[1]) }
        #
        return(root)
}

#
#MAIN
#

#
B0 = 10000
M  = 0.2
#
ff = .5*M #1 #0.5
#zeta = 45/((ff/M)+2) 
#gamma = 0.9#getGamma(zeta, ff, M)
n = 400
scope = 1
#tune scale on curvature
f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy
l = c()
for(i in seq(-scope, scope, length.out=30)){l=c(l, numDeriv::hessian(function(x){log(f(x))}, i))}#loglikelihood analogy
sd = median(sqrt(1/l), na.rm=T)
##h = hessian(z, 0-0.5, ff=ff, M=M)
#sd = sqrt(1/h)
##tune df on tails and thus emphasize or deemphasize the central peak.
#smaller df has wider tails and lower spread, higher df has thinner tails and more spread. 
gs = rt(n, 100)*sd  #rnorm(n, 0, sd) # #rnorm(n, 0, 2) #runif(n, -10, 10) #rt(n, 2)*sd+mu #
zs = z(gs, ff, M)
as = a(gs, ff, M)
bs = b(as, gs, M, B0)


##
#f = function(x){ ((ff/(M+ff))*(1-zeta)/(zeta*x)+1)^x - (M+ff)/M }
##uniroot(f, c(-1, -zeta-0.1))


#zz = getZeta(gamma, ff, M)
#alpha = getAlpha(gamma, ff, M)
#beta = getBeta(alpha, gamma, M, B0)
##print( PBar(alpha, beta, gamma, 0, M) )
#
##
#plot(SRR(0:B0, alpha, beta, gamma)-0:B0*M, type='l')
##curve(SRR(x, alpha, beta, gamma), from=0, to=B0)
#lines(0:B0, 0:B0*(ff))#+M))
#abline(v=zeta*B0)
