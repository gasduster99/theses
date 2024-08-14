rm(list=ls())

#
library(pracma)

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
#MAIN
#

#
cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(cpue)

#
#NOTE: N0Funk does not function correctly with $clone()
fit = prodModel$new(
        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #
        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
        gamma=2, #alpha=1, beta=cpue[1]/0.00049, gamma=2,   #parameters
        #lalpha=0, lbeta=log(cpue[1]/0.00049),    #reparameterize
        lalpha=-0.9395727, lbeta=7.8812027, 
	lq=log(0.00049), lsdo=-2.1176758 #log(0.01160256)  #nuisance parameters
        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
)

#K=2345.843, R=0.465
#pmLN = prodModel$new( dNdt=dNdt, time=1:TT, N0=3955.556, K=345.843, R=0.65, C=catch )
#pmLN$q = 0.00045
#pmLN$sdo  = 0.1 
#pmLN$model$observation = 'LN'
#optimize
optAns = fit$optimize(cpue, 
	c('lsdo', 'lalpha', 'lbeta'), 
	lower	= c(log(0.001), log(0.1), log(100)), 
	upper	= c(log(0.3), 0, log(10000)),  
	#gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
	#method 	= "Nelder-Mead", 
	cov	= T	
)

##dev.new()
#plot(cpue)
#fit$plotMean(add=T)
#fit$plotBand()

#
fitPT = prodModel$new(
        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #
        time=1:TT, catch=catch, #FtFmsy,    #schaefer   #constants
        alpha=1, beta=cpue[1]/0.00049, gamma=2,   #parameters
        lalpha=fit$lalpha, lbeta=fit$lbeta,    #reparameterize
        lq=fit$lq, lsdo=fit$lsdo  #nuisance parameters
        #xi=datGen$xi, zeta=datGen$zeta          #other incidentals to carry a
)

#optimize
optAns = fitPT$optimize(cpue,
        c('lsdo', 'lalpha', 'lbeta', 'gamma'),
        lower   = c(log(0.001), log(0.1), log(1000), eps()),
        upper   = c(log(0.3), 0, log(10000), 4),
        gaBoost = T, #list(maxiter=1e3, run=10, popSize=1e5), #T, #
        #method         = "Nelder-Mead", 
        cov     = T
)

##
#fitPT$plotMean(add=T, col="blue")
#fitPT$plotBand(col="blue")



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
