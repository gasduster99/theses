rm(list=ls())

#
library(lamW)
library(pracma)
library(mixtools)

#source('prodClass0.0.1.r')
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

##dNdt
#dNdt = function(t, y, rr, kk, catch){
#        #t      : current time step as requested by ode
#        #y      : value at previous time as requested by ode
#        #r      : growth rate given as numeric
#        #K/N0   : carrying capacity given as numeric
#	#C	: catch time series
#
#	#linearly interpolate catches
#        ft = floor(t)
#        q  = (t-ft)
#        Cl = catch[ft]
#        Cu = catch[ft+1]
#        C = q*Cu + (1-q)*Cl
#        if(q==0){ C=Cl } 
#
#        #
#	out = rr*y*(1-y/kk) - C[t]
#        R = exp(lalpha)*P*(1-exp(lbeta)*gamma*P)^(1/gamma)
#	out = R - C
#	#
#        return( list(out) )
#}

#
dPdt = function(t, P, lalpha, lbeta, gamma, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        #Fmsy = exp(lalpha)/gamma                                        #exp(lalpha)/
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta))^(gamma-1))        #exp(lalpha)*
        out = R - C #Fmsy*C*P
        #
        return( list(out) )
}

#
dPdtRPNat = function(t, P, Fmsy, lbeta, zeta, catch){ 
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        #Fmsy = exp(lalpha)/gamma 
	#zeta = boot::inv.logit(logitZeta
	gamma = getGamma(zeta)
	lalpha = log(getR(Fmsy, gamma))
	#
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta))^(gamma-1))        #exp(lalpha)*
        out = R - C #Fmsy*C*P
        #
        return( list(out) )
}

#
dPdtRP = function(t, P, lFmsy, lbeta, logitZeta, catch){ 
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        #Fmsy = exp(lalpha)/gamma 
	zeta = boot::inv.logit(logitZeta)
	gamma = getGamma(zeta)
	lalpha = log(getR(exp(lFmsy), gamma))
	#
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta))^(gamma-1))        #exp(lalpha)*
        out = R - C #Fmsy*C*P
        #
        return( list(out) )
}

#
getGamma = function(z){
	gammas = c(lambertW0(z*log(z))/log(z), lambertWm1(z*log(z))/log(z))
	return( gammas[1+(z>(1/exp(1)))] )
}

#
getZeta = function(g){
	(1/g)^(1/(g-1))
}

#
getR = function(Fmsy, g){
	Fmsy*g
}

#
getFmsy = function(r, g){
	r/g
}

##
#ptZeta = function(g){ (1/g)^(1/(g-1)) }
#ptXi = function(r, g){r/g}

#
#MAIN
#

#
cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(cpue)
#cpue  = c(0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
#catch = c(195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
#TT = length(cpue)

##
##NOTE: N0Funk does not function correctly with $clone()
#fit = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #
#        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
#        gamma=2, #alpha=1, beta=cpue[1]/0.00049, gamma=2,   #parameters
#        #lalpha=0, lbeta=log(cpue[1]/0.00049),    #reparameterize
#        lalpha=-0.9395727, lbeta=7.8812027, 
#	lq=log(0.00049), lsdo=-2.1176758 #log(0.01160256)  #nuisance parameters
#        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
#)
#
##K=2345.843, R=0.465
##pmLN = prodModel$new( dNdt=dNdt, time=1:TT, N0=3955.556, K=345.843, R=0.65, C=catch )
##pmLN$q = 0.00045
##pmLN$sdo  = 0.1 
##pmLN$model$observation = 'LN'
##optimize
#opt2 = fit$optimize(cpue, 
#	c('lsdo', 'lalpha', 'lbeta'), 
#	lower	= c(log(0.001), log(0.1), log(100)), 
#	upper	= c(log(0.3), 0, log(10000)),  
#	gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
#	#method 	= "Nelder-Mead", 
#	cov	= T,
#	fitQ = T
#)
#
##
#dev.new()
#plot(cpue)
#fit$plotMean(add=T)
#fit$plotBand()
#
##
#fitPT = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #
#        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
#        alpha=1, beta=cpue[1]/0.00049, gamma=0.8363341,   #parameters
#        lalpha=fit$lalpha, lbeta=fit$lbeta,    #reparameterize
#        lq=fit$lq, lsdo=fit$lsdo  #nuisance parameters
#        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
#)
#
##optimize
#optPT = fitPT$optimize(cpue,
#        c('lsdo', 'lalpha', 'lbeta', 'gamma'),
#        lower   = c(log(0.001), log(0.1), log(1000), eps()),
#        upper   = c(log(0.3), 0, log(10000), 4),
#        gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
#        #method         = "Neld er-Mead", 
#        cov     = T,
#	fitQ = T
#)
#
##
#fitPT$plotMean(add=T, col="blue")
#fitPT$plotBand(col="blue")
#
##
#fitPTRP = prodModel$new(
#        dNdt=dPdtRP, N0Funk=function(lbeta){exp(lbeta)}, #
#        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
#        #alpha=1, beta=cpue[1]/0.00049, gamma=0.8363341,   #parameters
#        #lalpha=fit$lalpha, lbeta=fit$lbeta,    #reparameterize
#        #dPdtRP = function(t, P, lFmsy, lbeta, logitZeta, catch)
#	lbeta=fitPT$lbeta, lFmsy=log(getFmsy(exp(fitPT$lalpha), fitPT$gamma)), logitZeta=logit(getZeta(fitPT$gamma)), 
#	lq=fit$lq, lsdo=fit$lsdo  #nuisance parameters
#        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
#)
#
##optimize
#optPTRP = fitPTRP$optimize(cpue,
#        c('lsdo', 'lFmsy', 'lbeta', 'logitZeta'),
#        lower   = c(log(0.001), log(0.1), log(1000),  logit(getZeta(0.1))),
#        upper   = c(log(0.3),   log(1), log(10000), logit(getZeta(4))),
#        gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
#        #method         = "Nelder-Mead", 
#        cov     = T,
#        fitQ = T
#)
##

#
fitPTRP = readRDS("PTRP.rda")

#
fitPTRPNat = prodModel$new(
        dNdt=dPdtRPNat, N0Funk=function(lbeta){exp(lbeta)}, #
        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
        #alpha=1, beta=cpue[1]/0.00049, gamma=0.8363341,   #parameters
        #lalpha=fit$lalpha, lbeta=fit$lbeta,    #reparameterize
        #dPdtRP = function(t, P, lFmsy, lbeta, logitZeta, catch)
       lbeta=fitPTRP$lbeta, Fmsy=exp(fitPTRP$lFmsy), zeta=boot::inv.logit(fitPTRP$logitZeta), 
       lq=fitPTRP$lq, lsdo=fitPTRP$lsdo  #nuisance parameters
        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
)

#
fitSRPNat = prodModel$new(
        dNdt=dPdtRPNat, N0Funk=function(lbeta){exp(lbeta)}, #
        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
        #alpha=1, beta=cpue[1]/0.00049, gamma=0.8363341,   #parameters
        #lalpha=fit$lalpha, lbeta=fit$lbeta,    #reparameterize
        #dPdtRP = function(t, P, lFmsy, lbeta, logitZeta, catch)
       lbeta=fitPTRP$lbeta, Fmsy=exp(fitPTRP$lFmsy), zeta=0.5, #boot::inv.logit(fitPTRP$logitZeta),
       lq=fitPTRP$lq, lsdo=fitPTRP$lsdo  #nuisance parameters
        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
)

#optimize
optSRPNat = fitSRPNat$optimize(cpue,
        c('lsdo', 'Fmsy', 'lbeta'),
        lower   = c(log(0.001), 0.1,  log(1000)),
        upper   = c(  log(0.3),   1, log(10000)),
        gaBoost = F, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
        #method         = "Nelder-Mead", 
        cov     = T,
        fitQ = T
)

##optimize
#optPTRPNat = fitPTRPNat$optimize(cpue,
#        c('lsdo', 'Fmsy', 'lbeta', 'zeta'),
#        lower   = c(log(0.001), 0.1,  log(1000), getZeta(0.1)),
#        upper   = c(  log(0.3),   1, log(10000), getZeta(4)  ),
#        gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
#        #method         = "Nelder-Mead", 
#        cov     = T,
#        fitQ = T
#)

#1, 4
set.seed(4)

##
#samNat = rmvnorm(1000, c(fitPTRPNat$Fmsy, fitPTRPNat$zeta), fitPTRPNat$rsCov[c("Fmsy", "zeta"), c("Fmsy", "zeta")])
#png("namibibRP.png")
##samE = cbind( exp(sam[,1]), boot::inv.logit(sam[,2]) )
#plot(samNat[,1], samNat[,2], xlab="Fmsy", ylab="Bmsy/B0", main="Namibian Hake RP Estimates", xlim=c(0.1, 0.7), ylim=c(0.1, 0.8), col="white")
##point(c(fitPTRP$lFmsy, fitPTRP$logitZeta), pch=19)
#ellipse(colMeans(samNat), cov(samNat), alpha=0.5)
#points(t(colMeans(samNat)), pch=19)
#sam2 = rnorm(1000, fit$lalpha, fit$rsCov[c("lalpha"), c("lalpha")])
#points(exp(fit$lalpha)/2, 0.5, pch=19, col='red')
#segments(quantile(exp(sam2)/2, 0.25), 0.5, quantile(exp(sam2)/2, 0.75), 0.5, col='red')
#abline(h=0.5)
#legend("topright", legend=c("Three-Parameter PT RPs", "Two-Parameter Schaefer RPs"), pch=19, col=c('black', 'red'))
#dev.off()

#
fitPTRPNat = readRDS("PTRPNat.rda")
samNat = rmvnorm(1000, c(fitPTRPNat$Fmsy, fitPTRPNat$zeta), fitPTRPNat$rsCov[c("Fmsy", "zeta"), c("Fmsy", "zeta")])
#png("namibibRP.png")
#samE = cbind( exp(sam[,1]), boot::inv.logit(sam[,2]) )
plot(samNat[,1], samNat[,2], xlab="Fmsy", ylab="Bmsy/B0", main="Namibian Hake RP Estimates", xlim=c(0.1, 0.7), ylim=c(0.1, 0.8), col="white")
#point(c(fitPTRP$lFmsy, fitPTRP$logitZeta), pch=19)
ellipse(colMeans(samNat), cov(samNat), alpha=0.5)
points(t(colMeans(samNat)), pch=19)
sam2 = rnorm(1000, fit$lalpha, fit$rsCov[c("lalpha"), c("lalpha")])
points(exp(fit$lalpha)/2, 0.5, pch=19, col='red')
segments(quantile(exp(sam2)/2, 0.25), 0.5, quantile(exp(sam2)/2, 0.75), 0.5, col='red')
abline(h=0.5)
legend("topright", legend=c("Three-Parameter PT RPs", "Two-Parameter Schaefer RPs"), pch=19, col=c('black', 'red'))
#dev.off()

#
#RP Plot
#

#
fitPTRP = readRDS("PTRP.rda")
fitPT = readRDS("PT.rda")
fit = readRDS("schaefer.rda")

##1, 4
#set.seed(4)
##
#sam = rmvnorm(100, c(fitPTRP$lFmsy, fitPTRP$logitZeta), fitPTRP$rsCov[c("lFmsy", "logitZeta"), c("lFmsy", "logitZeta")])
#samE = cbind( exp(sam[,1]), boot::inv.logit(sam[,2]) )
#plot(exp(sam[,1]), boot::inv.logit(sam[,2]), col='white', xlab="Fmsy", ylab="B*/B0", main="Namibian Hake RP Estimates")
##point(c(fitPTRP$lFmsy, fitPTRP$logitZeta), pch=19)
#ellipse(colMeans(samE), cov(samE), alpha=0.5)
#points(t(colMeans(samE)), pch=19)
#sam2 = rnorm(1000, fit$lalpha, fit$rsCov[c("lalpha"), c("lalpha")])
#points(exp(fit$lalpha)/2, 0.5, pch=19, col='red')
#segments(quantile(exp(sam2)/2, 0.25), 0.5, quantile(exp(sam2)/2, 0.75), 0.5, col='red')
#abline(h=0.5)

##
#plot(ptXi(exp(fitPT$lalpha),fitPT$gamma), ptZeta(fitPT$gamma), xlim=c(0, 3), ylim=c(0,1))
#abline(h=0.5)



##
#dev.new()
#sam = pmLN$plotRS(m=1e5)#save='test.jpg')

##
#pmN = prodModel$new( dNdt=dNdt, time=1:TT, N0=3955.556, K=3955.556, R=0.38, C=catch )
#pmN$q = 0.00045
#pmN$sdo = 0.1
#pmN$model$observation = 'N'
##
#optAns = pmN$optimize(cpue, 
#	c('sdo', 'R', 'K'), 
#	lower	= c(0.001, 0, 0), 
#	upper	= c(0.3, 1, 1e5),  
#	gaBoost = T, 
#	cov	= T	
#)
#
##
#pmN$plot(add=T, col='red')
