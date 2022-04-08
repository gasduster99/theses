rm(list=ls())

#
library(GA)
library(numDeriv)
#library(rootSolve)

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

#number of samples in design
n = 50 #500 maxed out method

#LHS boundaries
xiMin = 0
xiMax = 4 #2 #3.5
zetaMin = 0.1 #0.1
zetaMax = 0.7 #0.3 #0.7 #1
#zeta range should not be too small

#xi and zeta Bins (bin defined by left edge right egde not used) [n+1 used to give n left edges]
xiL = seq(xiMin, xiMax, length.out=n+1)
zetaL = seq(zetaMin, zetaMax, length.out=n+1)
#
Llist = rep(T, n) #a list to track which zeta bin still needs a sample
zetaList = rep(NA, n) #a list the sampled zeta value in each bin
xiList = rep(NA, n)   #a list the sampled xi value in each bin

#details of how to sample the curvy bit for comuting sd
scope = 2 #how far to sample the inverse function for curvature
scopeN = 30 #how many samples to take of the inverse function curvature
scopeRun = seq(-scope, scope, length.out=scopeN)
#start at fudge (hessian sd) slowing increase fudge by fudgeFactor so that we eventually search flat region, but still allow for repetition at calculated sd 
fudge = 1
giveUp = 1000
fudgeFactor = 0.01
#start df at dfStart with thin tails (to explore curvature), decrement by dfFactor to transition to explore flat regiem 
dfStart = 100
df = dfStart
dfFactor = 1 #dfFactor, dfStart, and fudgeFactor are balanced to first decrease df, and then crank on sd

#
xii = 1
xiShuffle = sample(xiL[-(n+1)])
while( any(Llist) ){
	#for each xi column (n of them)
	xi = runif(1, xiShuffle[xii], xiL[which(xiL==xiShuffle[xii])+1])
	ff = xi*M 
	
	#tune scale on curvature
	f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy	
	hesses = sapply(scopeRun, function(i){numDeriv::hessian(function(x){log(f(x))}, i)})	
	sd = median(sqrt(1/hesses), na.rm=T)
	#may (maybe not) be beneficial to consider adding a tuning step for the flatter sections of getZeta(z)??

	##tune df on student t tails and thus emphasize or de-emphasize the central peak.
	#smaller df has wider tails and lower spread (since it samples flat regiem of getZeta), 
	#higher df has thinner tails and more spread (since is samples the spiky bit of getZeta). 
	gs = rt(1, df)*sd*fudge #rt(1, 100)*sd*fudge  
	
	#propose sample of zeta
	zs = z(gs, ff, M)
	#reject: catch the case where zs is propsed outside of z range
	if(zs<zetaMin | zs>zetaMax){ next } 
	
	#NOTE: compute auxilliary quantities
	#NOTE: make sure that reverse a,b,g->xi,zeta and xi,zeta->a,b,g calcs make sense; else reject	
	#as = a(gs, ff, M)
	#bs = b(as, gs, M, B0)
	#if( reverse a,b,g->xi,zeta and xi,zeta->a,b,g fail ){ next } #reject

	#find left edge of bin sampled
	bin = max(which(zetaL<=zs))
	#if sampled bin not yet sampled, then save
	if( Llist[bin] ){
		#record values
		xiList[bin] = xi 
		zetaList[bin] = zs
		#
		#NOTE: generate data and catelogue design elements in a model object in the design directory
		#NOTE: catelogue by bin lables (left edges of xiL, zetaL)s
		#
		#knobs for making the design
		Llist[bin] = F 	#turn off bin
		xii = xii+1 	#increment xi counter
		df = 100    	#reset df to thin tails
		fudge = 1   	#reset fudge to thin tails
	}else{ # catch case where hyperbola gap requires heavy tails
		df = max(1, df-1)  #decrement df to broadend search to get the center bald spot
		fudge = fudge+0.01 #increase fudge factor to get a wider sd
		if( fudge>giveUp ){ break } #give up (it mostly just misses one of two bins so its usually fine)
	}
}

#
png('design.png')
plot(xiList, zetaList, pch=20)
abline(v=xiL, lty=3)
abline(h=zetaL, lty=3)
dev.off()

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
