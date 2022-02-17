rm(list=ls())

#
library(mvtnorm)
library(bayesmeta)
library(matrixStats)
library(foreach)
library(doParallel)

#
#FUNCTIONS
#

#
pDiff = function(t1, t2){ (t1-t2)*(t1>t2) }
makeXX = function(){ cbind(1, sapply(time[-TT], FUN=pDiff, t1=time))}

#
ff = function(t){
        x = cbind(1, t(sapply(time[-TT], FUN=pDiff, t1=t)))
        x%*%bStar
}
ff = Vectorize(ff)

#
makeXY = function(ts, knots, asMat=T){
        #       
        x = c()
        for(t in ts){
                tauStar = ceiling(t)-1 #t-1 #floor(t)
                x = rbind( x, c( 1, (t^2/2-t*knots)*(t>knots) - (tauStar^2/2-tauStar*knots)*(tauStar>=knots)) )
                #B_0: 1 -or- t-tauStar
        }
	if( asMat ){ return(x) }
        return(as.data.frame(x))
}

#
likelihood = function(y, X, beta, obsVar, log=F){
        #y: an n vector of responses
        #X: an nxp matrix of predictors
        #beta: a p vector of parameters
        #log: should the function return the density or log density 
        #
        #value: evaluation of the likelihood (or log likelihood)

        #
	like = dmvnorm(y, X%*%beta, obsVar, log=log)
        #like = dbinom(y, m, piFunk(X%*%beta), log=log)
        out = sum(like)
        if( !log ){ out=prod(like) }
        #
        return(out)
}

#
betaPrior = function(beta, v, log=F){
	#beta: a vector of the slopes
	#v: a scalar of the pooled variance among the betas
	#
	#value: evaluation of the beta prior (or its log)
	
	#
	bLike = dnorm(beta, sqrt(v), log=log)
	out = sum(bLike)
	if( !log ){ out=prod(bLike) }
	#
	return(out)
}

#
vPrior = function(v, log=F){
        #v: a scalar of the pooled variance among the betas
        #
        #value: evaluation of the v prior (or its log)

        # 
        vP = dhalfcauchy(v, 1, log=log)
        #
        return(out)
}


#
#DATA
#

#                                           606
catch = c(94, 212, 195, 383, 320, 402, 366, 505, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(catch)
time = 1:TT

#
X = makeXY(time, time[-TT])
XX = makeXX()

#
#OPTIMIZE
#

#lmFit = lm(catch~X-1, wieght=1/catch^2)
#pCov = summary(lmFit)$cov.scaled

#snr=1/cv=m/sd
snr = 100000 #2 #3 #(0.1, 100000)
sig = diag(catch)/snr #100000
ISig = solve(sig)
coefCov = solve(t(X)%*%ISig%*%X)
coefWLS = coefCov %*% t(X)%*%ISig%*%catch
#
lvSD = 1
covP = 0.8 #1.5

#
vScale = 1#snr

#
#SAMPLE
#

#
M = 1000000

#
threads=8
registerDoParallel(cores=threads)
post = foreach(c=1:threads)%dopar%{
	#
	ff = function(t){
	        x = cbind(1, t(sapply(time[-TT], FUN=pDiff, t1=t)))
	        x%*%bStar
	}
	ff = Vectorize(ff)

	#
	intPred = matrix(NA, nrow=M+1, ncol=TT); intPred[1,]=XX%*%coefWLS
	bluePred = matrix(NA, nrow=M+1, ncol=TT); bluePred[1,]=X%*%coefWLS
	catchPred = matrix(NA, nrow=M+1, ncol=TT); catchPred[1,]=XX%*%coefWLS
	#
	beta0Post = matrix(NA, nrow=M+1, ncol=1); beta0Post[1]=coefWLS[1]
	betaPost = matrix(NA, nrow=M+1, ncol=TT-1); betaPost[1,]=coefWLS[-1]
	vPost = matrix(NA, nrow=M+1, ncol=1); vPost[1]=1

	#
	for(m in 1:M){
		#PROPOSE
		bPost = c(beta0Post[m], betaPost[m,])	
		bStar = t(rmvnorm(1, bPost, coefCov*covP)) 
		vStar = exp(rnorm(1, log(vPost[m]), lvSD))

		#
		top = likelihood(catch, X, bStar, sig, log=T) + 
			betaPrior(bStar[-1], vStar, log=T)    + 
			dhalfcauchy(vStar, vScale, log=T)     + 
			log(vStar)
		bot = likelihood(catch, X, bPost, sig, log=T) + 
			betaPrior(bPost[-1], vPost[m], log=T) + 
			dhalfcauchy(vPost[m], vScale, log=T)  + 
			log(vPost[m])
		q = min(1, exp(top-bot))
	
		#reject
	        beta0Post[m+1,] = beta0Post[m,]
	        betaPost[m+1,]  = betaPost[m,]
		vPost[m+1,]     = vPost[m]
		#
		catchPred[m+1,] = catchPred[m,]
	        bluePred[m+1,]  = bluePred[m,] #rmvnorm(1, bluePred[m,], sig)
		intPred[m+1,]   = intPred[m,]
		#accept
	        if( rbinom(1, 1, q) ){
			#	
			beta0Post[m+1,] = bStar[1]
			betaPost[m+1,]  = bStar[-1]
			vPost[m+1,]     = vStar
			#
			catchPred[m+1,] = XX%*%bStar
			bluePred[m+1,]  = rmvnorm(1, X%*%bStar, sig)
	        	intPred[m+1,]   = sapply(time, FUN=function(x){integrate(ff, x-1, x)$value})
		}

		##thin & save    
                #if( !m%%thinEvery ){
                #        #thinned index
                #        mm = m/thinEvery

                #        #psi
                #        psi[mm,] = psiNow
                #        #t2|psi
                #        t2[mm]   = rinvgamma(1, (n-k)/2+at-1, fbNow$S2/2+bt)
                #        #beta|t2, psi
                #        beta[mm,]= rmvnorm(1, fbNow$bHat, t2[mm]*strongInv(t(fbNow$LIX)%*%fbNow
                #}
	}
        #
	flame = seq(1000, M+1, 1)
        return( list(beta0Post=beta0Post[flame], betaPost=betaPost[flame,], vPost=vPost[flame], catchPred=catchPred[flame,], bluePred=bluePred[flame,], intPred=intPred[flame,]) )
}

#
beta0Post = sapply(post, function(p){p['beta0Post']})
beta0Post = do.call(rbind, beta0Post)
betaPost = sapply(post, function(p){p['betaPost']})
betaPost = do.call(rbind, betaPost)
vPost = sapply(post, function(p){p['vPost']})
vPost = do.call(rbind, vPost)
catchPred = sapply(post, function(p){p['catchPred']})
catchPred = do.call(rbind, catchPred)
bluePred = sapply(post, function(p){p['bluePred']})
bluePred = do.call(rbind, bluePred)
intPred = sapply(post, function(p){p['intPred']})
intPred = do.call(rbind, intPred)


#burn = 100
thin = 1000
thinner = seq(1, length(beta0Post), thin)
#
beta0Post = beta0Post[thinner]
betaPost  = betaPost[thinner,]
vPost     = vPost[thinner]
#
catchPred = catchPred[thinner,]
bluePred = bluePred[thinner,]

#
betaHat = c(mean(beta0Post), colMeans(betaPost))
#XX = makeXX()

##
#png('ridgeSpline.png')
#boxplot(intPred, outline=F, border='red')
#points(time, catch, pch=20)
#i = 1
#for(sd in sqrt(diag(sig))){ segments(time[i], catch-1.96*sd, time[i], catch+1.96*sd); i=i+1 }
#lines(time, colMeans(catchPred), lwd=3)
#polygon(c(time, rev(time)), c(colQuantiles(catchPred, probs=0.025), rev(colQuantiles(catchPred, probs=0.975))),
#	col=adjustcolor('black', alpha.f=0.4),
#	border=NA
#)
#lines(time, colMeans(bluePred), lwd=3, col='blue')
#dev.off()
#
##lines(time, colQuantiles(catchPred, probs=0.975))
##lines(time, colQuantiles(bluePred, probs=0.025), col='blue')
##lines(time, colQuantiles(bluePred, probs=0.975), col='blue')

#
ff = function(t){
        x = cbind(1, t(sapply(time[-TT], FUN=pDiff, t1=t)))
        x%*%betaHat
}
ff = Vectorize(ff)
#points(time, sapply(time, FUN=function(x){integrate(ff, x-1, x)$value}), col='red')

#
X2 = sapply(time[-TT], FUN=pDiff, t1=time)
X4 = matrix(0, nrow=TT, ncol=TT)
for(i in time){
        X4[i, 1] = 1 #time[i]-(time[i]-1)
        if(i-1<1){ next }
        for(j in 1:(i-1)){
                X4[i, j+1] = ((time[i]^2)/2-time[i]*time[j]) - ((time[i-1]^2)/2-time[i-1]*time[j])
        }
}
#
sp4 = lm(catch~X4-1)
coef = sp4$coef

#
#INTERPOLATION PLOTS
#

#
png('intCurvesSmallCV.png')
#
plot(time, cbind(1, X2)%*%sp4$coefficients, type='l', lwd=2, xlab="Time", ylab="Catch", ylim=c(0, 1300), col='blue', main="Interpolated Instantaneous Catch")
#boxplot(intPred, outline=F, border='red')
#points(time, catch, pch=20)
#i = 1
#for(sd in sqrt(diag(sig))){ segments(time[i], catch-1.96*sd, time[i], catch+1.96*sd); i=i+1 }
lines(time, colMeans(catchPred), lwd=3)
polygon(c(time, rev(time)), c(colQuantiles(catchPred, probs=0.025), rev(colQuantiles(catchPred, probs=0.975))),
	col=adjustcolor('black', alpha.f=0.4),
	border=NA
)
legend('topleft', legend=c('Observations Assumed Known', 'Observations Assumed Uncertain'), col=c('blue', 'black'), lwd=c(2, 3))
#lines(time, colMeans(bluePred), lwd=3, col='blue')
dev.off()

#
#DATA PLOTS
#

#
png('dataPlotsSmallCV.png')
#
plot(time, X4%*%coef, 'l', col='blue', ylim=c(0, 700), lwd=2, ylab="Catch", xlab="Time", main="Observed Catch with Predictive Interpolations") #, ylim=c(0, 1500))
points(time, catch, pch='-', cex=1.75)#20)

#
i = 1
for(sd in sqrt(diag(sig))){ segments(time[i], catch[i]-1.96*sd, time[i], catch[i]+1.96*sd); i=i+1 }

##
#lines(time, colMeans(catchPred), lwd=3)
#polygon(c(time, rev(time)), c(colQuantiles(catchPred, probs=0.025), rev(colQuantiles(catchPred, probs=0.975))),
#	col=adjustcolor('black', alpha.f=0.4),
#	border=NA
#)
#lines(time, colMeans(bluePred), lwd=3, col='blue')

#
polygon(c(time, rev(time)), c(colQuantiles(bluePred, probs=0.025), rev(colQuantiles(bluePred, probs=0.975))),
        col=adjustcolor('black', alpha.f=0.4),
        border=NA
)
lines(time, colMeans(bluePred), lwd=3, col='black')
legend('topleft', legend=c('Observations Assumed Known', 'Observations Assumed Uncertain'), col=c('blue', 'black'), lwd=c(2, 3))
dev.off()









#
#JUNK
#

##
		#top = likelihood(Y, mu, log=T) + muPrior(mu, lambdaStar, gammaStar, log=T) + lambdaPrior(lambdaStar, pp, log=T) + log(lambdaStar)
	        #bot = likelihood(Y, mu, log=T) + muPrior(mu, lambda[m], gammaM, log=T) + lambdaPrior(lambda[m], pp, log=T) + log(lambda[m])
	        #q = min(1, exp(top-bot))
	
	        ##reject
	        #beta[m+1,] = beta[m,]
	        #lambda[m+1] = lambda[m]
	        ##accept
	        #if( rbinom(1, 1, q) ){
	        #        beta[m+1,] = betaStar
	        #        lambda[m+1] = lambdaStar
	        #}

