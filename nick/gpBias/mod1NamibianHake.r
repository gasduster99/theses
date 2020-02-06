rm(list=ls())

#
library(geoR)
library(boot)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, Cs, lhs, gamma, H, delta){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #lhs  	: logit of h* given as numeric
        #Cs   	: C* given as numeric
        #gamma  : recruitment parameter given as numeric

        #Derivative
        C = catch[t]
        #
	hs = inv.logit(lhs)
        sh = (1-delta)*(1-hs)
        alpha = (1-sh)/(1-hs) * (1+(gamma*hs)/(1-sh))^(1/gamma)
        beta = (hs^2)/((1-hs)*(1-sh+gamma*hs)*Cs)
        #
        R = alpha*(y-C)*(1-beta*gamma*(y-C))^(1/gamma)
        out = R - delta*(y-C) - C
        #
        return( list(out) )
}

#
N0Funk = function(delta, Cs, lhs, gamma){
	#parameterization
	hs = inv.logit(lhs)
	sh = (1-delta)*(1-hs)
        alpha = (1-sh)/(1-hs) * (1+(gamma*hs)/(1-sh))^(1/gamma)
        beta = (hs^2)/((1-hs)*(1-sh+gamma*hs)*Cs)
	#Virgin 
	r0 = delta
	R0 = (1-(r0/alpha)^gamma)*r0/(beta*gamma)
	P0 = R0/r0
	#
	return( P0 )
}

#
#FIT
#

#
cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(cpue)
#
M = 0.2

##initialize
#mod = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=1-exp(-M), catch=catch, 
#	Cs    = 263.9692, 	   #mean(catch),
#	lhs   = logit(0.1702508),
#	gamma = 1,  
#	lq    = log(0.0004271659), #q = 0.00047, 
#	sdo   = 0.146974
#)
#mod$model$observation = 'LN'

#
mod = readRDS('modAll.rds')
mod$gamma = .Machine$double.eps
mod$lq    = log(1e-5)
mod$lhs   = logit(0.01)
mod$Cs    = 470.0084 

#
dev.new()
plot(cpue, ylim=c(0, max(cpue)*1.1))
mod$iterate()
mod$plotMean(add=T, col='red')
#optimize
optAns = mod$optimize(cpue,
        c('sdo', 'Cs', 'lhs', 'lq'),
        lower   = c(0.001, min(catch)*0.5, logit(0.001), log(1e-7)),
        upper   = c(1, max(catch), logit(0.5), log(1e-2)),
        gaBoost = list('parallel'=8), #T,
        cov     = T
)
mod$plotMean(add=T)
mod$plotBand()
#
mod$printSelf()
mod$save('gammaZero.rds')


#
#TRANSFORM
#

#
boxTrans = function(x, lambda, lambda2=0){
        if(lambda==0){return(log(xi+lambda2))}
        return( ((x+lambda2)^lambda-1)/lambda )
}

#hats
hsHat = inv.logit(mod$lhs)
FsHat = -log(1-hsHat)
xiHat = FsHat/M
zetaHat = 1/(xiHat+2)

#rs distributions
dev.new()
sam  = mod$plotRS()
hs   = inv.logit(sam[,'lhs'])
Fs   = -log(1-hs)
xi   = Fs/M
zeta = 1/(xi+2)
gamma = logit(2*zeta)

#ins = gamma-min(gamma)+.Machine$double.eps
out = boxcoxfit(gamma, lambda2=T)


















##
##dNdt
#dNdtAB = function(t, y, alpha, beta, gamma, catch, delta){
#        #t      : current time step as requested by ode
#        #y      : value at previous time as requested by ode
#        #alpha  : recruitment parameter given as numeric
#        #beta   : recruitment parameter given as numeric
#        #gamma  : recruitment parameter given as numeric
#
#	#Derivative
#	C = catch[t]
#	R = alpha*(y-C)*(1-beta*gamma*(y-C))^(1/gamma)
#	out = R - delta*(y-C) - C
#	#
#	return( list(out) )
#}

##initialize
#mod = prodModel$new(dNdt=dNdt, time=1:TT, N0Funk=N0Funk, delta=1-exp(-M), catch=catch, 
#	alpha = 1.99573,
#	beta  = 5.904659e-06,
#	gamma = 1, 
#	q     = 0.0007027444,  
#	sdo   = 0.3548753
#)
#mod$model$observation = 'N'
#c('sdo', 'alpha', 'beta', 'q')
#c(0.001, mod$delta, .Machine$double.eps, 1e-5)
#c(1, 10, 1, 1e-2)







#datGen$plotMean()
##tuning parameter
#kappa = 0
#w1    = 0.001
#winf  = 0.001







##Virgin
#R = (1-(rh/alpha)^gamma)*rh/(beta*gamma)
#
##Recruitment
#if(t>k){
#	R = alpha*(y-C[t-1])*(1-beta*gamma*(y-C[t-1]))^(1/gamma) 
#}
#out = 
#C = h*P
#out = R*y*(1-y/K) - C[t]
##
#return( list(out) )
