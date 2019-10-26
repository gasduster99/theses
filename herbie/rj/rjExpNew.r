rm(list=ls())

library(ggplot2)

load('rjMCMC10000.RData') #load('rjMCMC100000.RData')#

kNum = Mode(k)#2 #5

#DATA
hist(classAVG,
	breaks=7,
	xlab='Score',
	ylab='',
	main='Data Histogram',
	freq=F
)

#TRACE PLOTS
dev.new(width=15, height=7)
layout(matrix(c(1,2,3), nrow=1, ncol=3))
#weights
matplot(wO[[kNum]], 
	type='l',
	xlab='Sample',
	ylab=expression('w'[j]),
	main=expression('w')
)

#means
matplot(muO[[kNum]], 
	type='l',
	xlab='Sample', 
	ylab=expression(mu[j]),
	main=expression(mu)
)

#variances
matplot(sig2O[[kNum]], 
	type='l', 
	xlab='Sample', 
	ylab=expression(sigma[j]^2),
	main=expression(sigma^2),
	ylim=c(0,1)
)

dev.new(width=10, height=7)
layout(matrix(c(1,2), nrow=1, ncol=2))

#components
plot(	k,	
	type='l',
	xlab='Sample',
	ylab='k',
	main='k'
)


#beta
plot(	beta, 
	type='l',
	xlab='Sample',
	ylab=expression(beta),
	main=expression(beta)
)


#JOINT ESTIMATION
dev.new(width=7, height=5)#10)
#par(mfrow=c(2,1), mar=c(5, 4, 4, 15))#12))

#likelyhood
lMaxi = which(likelyHood==max(likelyHood))
lW = wRJ[[lMaxi]]
lMu = muRJ[[lMaxi]]
lSig2 = sig2RJ[[lMaxi]]
lK = k[lMaxi]
lBeta = beta[lMaxi]

plot(likelyHood,
	type='l',
	xlab='Sample #',
	ylab='',
	main=expression('Evaluations of the Log-Likelyhood')
)
points(lMaxi, likelyHood[lMaxi], col='red')
mtext(bquote(atop(atop(atop(atop(textstyle(hat(k)==.(lK)),
		textstyle(hat(w)%~~%bgroup('[',.(toString(round(lW, 2))),']'))),
		textstyle(hat(mu)%~~%bgroup('[',.(toString(round(lMu, 2))),']'))),
		textstyle(hat(sigma^2)%~~%bgroup('[',.(toString(round(lSig2, 4))),']'))),
		textstyle(hat(beta)%~~%.(toString(round(lBeta, 4)))))			
	),
	side=4,
	line=3,
	las=2,
	cex.lab=1
)




##posterior
#pMaxi = which(propPost==max(propPost))
#pW = wRJ[[pMaxi]]
#pMu = muRJ[[pMaxi]]
#pSig2 = sig2RJ[[pMaxi]]
#pK = k[pMaxi]
#pBeta = beta[pMaxi]
#
#dev.new()
#par(mar=c(5, 4, 4, 15))
#plot(propPost,
#	type='l',
#	xlab='Sample #',
#	ylab='',
#	main=expression('Evaluations of '*p(k,w,z,mu,sigma^2))
#)
#points(pMaxi, propPost[pMaxi], col='red')
#mtext(bquote(atop(atop(atop(atop(textstyle(hat(k)==.(pK)),
#		textstyle(hat(w)%~~%bgroup('[',.(toString(round(pW, 2))),']'))),
#		textstyle(hat(mu)%~~%bgroup('[',.(toString(round(pMu, 2))),']'))),
#		textstyle(hat(sigma^2)%~~%bgroup('[',.(toString(round(pSig2, 4))),']'))),
#		textstyle(hat(beta)%~~%.(toString(round(pBeta, 4)))))			
#	),
#	side=4,
#	line=3,
#	las=2,
#	cex.lab=1
#)



#DISTRIBUTIONS
WO = data.frame(value=numeric(0))
MUO = data.frame(value=numeric(0))
SIG2O = data.frame(value=numeric(0))
for (j in seq(1, kNum)){
	dfW = data.frame(value=wO[[kNum]][,j]) 
	dfW$Component = sprintf('%d',j)
	WO = rbind(WO, dfW)

	dfMu = data.frame(value=muO[[kNum]][,j]) 
	dfMu$Component = sprintf('%d',j)
	MUO = rbind(MUO, dfMu)

	dfSig2 = data.frame(value=sig2O[[kNum]][,j]) 
	dfSig2$Component = sprintf('%d',j)
	SIG2O = rbind(SIG2O, dfSig2)
}

#weights
dev.new()
print( ggplot(WO, aes(value, fill=Component)) + geom_histogram(alpha=0.6, position='identity') + opts(title=expression('w')) )


#means
dev.new()
print( ggplot(MUO, aes(value, fill=Component)) + geom_histogram(alpha=0.6, position='identity') + opts(title=expression(mu)) )


#variances
dev.new()
print( ggplot(SIG2O, aes(value, fill=Component)) + geom_histogram(alpha=0.6, position='identity') + opts(title=expression(sigma^2)) )


#components
dev.new()
hist(	k,
	xlab='k',
	ylab='',
	main='k'	
)


##beta
#dev.new()
#hist(	beta,
#	xlab=expression(beta),
#	ylab='',
#	main=expression(beta)	
#)














