rm(list=ls())

#
library(pracma)
#
source('prodClass0.1.1.r')
source('gpClass0.0.1.r')

#
#FUNCTIONS
#

#
dNdt = function(t, y, rr, K, R){ list(rr*y*(1-y/K) - R[t]) }

#
S2 = function(X0, X1, s2, v){
        maxD = dnorm(X0, X0, sqrt(v), log=T)
        s2*mapply(function(x0, m){
                exp(dnorm(X1, x0, sqrt(v), log=T)-m)
        }, X0, maxD) 
}


#
#SOUTH KOREA
#

#
skN = 51.47*10^6
skDat = read.table('sk.csv', header=T, sep=',')
skR = rowSums(skDat[,3:4])
skT = nrow(skDat)

#
sk = prodModel$new(time=1:skT, dNdt=dNdt, R=skR, model=list(observation="N"),
	rr   = 0.3994632, 
	K    = 7523.6, 
	lsdo = 4.498353, #-1,#
	N0   = 19.998  
)

#
skOut = sk$optimize(skDat$activeCases, 
	c('lsdo', 'rr', 'K', 'N0'),
	lower   = c(0, 0.3, 10^3, 5),
	upper   = c(10, 1, 10^5, 50),       
	#gaBoost = list(run=30, popSize=1e4, parallel=8),
        cov     = T
)
#
writeLines('South Korea')
sk$printSelf()


#
#ITALY
#

#
itN = 60.48*10^6
itDat = read.table('italy.csv', header=T, sep=',')
itR = rowSums(itDat[,3:4])
itT = nrow(itDat)

#
it = prodModel$new(time=1:itT, dNdt=dNdt, R=itR, model=list(observation="N"),
        rr   = 0.2713875, #0.3179662,
        K    = 42572.68,  #20129.5,
        lsdo = 5.277671, #5.172899, #5.216676,
        N0   = 27.83871 #13.80569
)

#
itOut = it$optimize(itDat$activeCases, 
	c('lsdo', 'rr', 'K', 'N0'),
	lower   = c(0, 0.26, 10^3, 5),
	upper   = c(10, 1, 10^5, 50),       
	#gaBoost = list(run=30, popSize=1e4, parallel=8),
        cov     = T
)
#
writeLines('\nItaly')
it$printSelf()

#
#USA
#

#-127.3914
usN = 60.48*10^6
usDat = read.table('us.csv', header=T, sep=',')
usR = rowSums(usDat[,3:4])
usT = nrow(usDat)

#
us = prodModel$new(time=1:usT, dNdt=dNdt, R=usR, model=list(observation="N"),
        rr   = 0.2747354,	#0.2901101,	
        K    = 777867.2,	#232874.5, 
        lsdo = 3.464444,	#3.026116, 	#4.866844, 	#5.216676,
        N0   = 2.224132		#1.686993		#12		#13.80569
)

#
usOut = us$optimize(usDat$activeCases, 
	c('lsdo', 'rr', 'K', 'N0'),
	lower   = c(0, 0.25, 10^4, 1),
	upper   = c(10, 1, 10^6, 50),       
	gaBoost = list(run=30, popSize=1e4, parallel=8),
        cov     = T
)
#
writeLines('\nUSA')
us$printSelf()

##
#gp = gpModel$new( S2=S2, X=us$time, Y=us$R, 
#	v  = 1,
#	s2 = var(us$R),
#	B  = as.vector(lm(us$R~us$time)$coef),
#	g  = 0#,
#	#obsV = c(rep(20, 18), rep(100, 12))
#)
#bnd = rep(Inf, length(gp$B))
#optOut = gp$fit(c('v', 's2', 'g'),#, 'B'),
#        lower   = c(eps(), eps(), -Inf),#, -bnd),
#        upper   = c(Inf, Inf, Inf),#, bnd),
#        cov     = T
#)
##
#tXX = seq(min(gp$X), max(gp$X)+15, 0.1) #seq(max(gp$X)+1, max(gp$X)+60, 1)
#pm = gp$predictMean(tXX)
#psd = sqrt(diag(gp$predictVar(tXX)))
#
##
#plot(gp$X, gp$Y, xlim=c(min(tXX), max(tXX)), ylim=c(0, max(pm)))
#lines(tXX, pm)
#lines(tXX, qnorm(0.025, pm, psd))
#lines(tXX, qnorm(0.975, pm, psd))

##
#plot(usDat$activeCases, col='black', pch=19, ylim=c(0, 9*10^5), xlim=c(0, max(tXX)))
#us$plotMean(add=T, col='black')
#us$plotBand()
##
#us$time = c(us$time, tXX)
#us$R = c(us$R, pm)
#us$iterate()
#lines(tXX-1, us$N[tXX-1])
#abline(h=us$K)
#abline(h=qnorm(0.975, us$K, sqrt(us$rsCov['K', 'K'])), lty=3)
#abline(h=qnorm(0.025, us$K, sqrt(us$rsCov['K', 'K'])), lty=3)


##
#plot(us$R, xlim=c(min(gp$X), max(tXX)) )
#lines(tXX, pm)
#lines(tXX, pm+1.96*psd, lty=2)
#lines(tXX, pm-1.96*psd, lty=2)

#
#PLOT
#

#
jpeg('cases.jpg')
#
plot(itDat$activeCases, 
	main = 'Corona Virus Cases Thru Time',
	ylab = '# Active Cases',
	xlab = 'Days Since Feb 14th 2020',
	ylim = c(0, max(c(skDat$activeCases, itDat$activeCases, usDat$activeCases))),
	col  = 'darkred', 
	pch  = 19
)
it$plotMean(add=T, col='red')
it$plotBand(col='red')
#abline(h=it$K, col='red')
#abline(h=qnorm(0.975, it$K, it$rsCov['K', 'K']), col='red', lty=2)
#abline(h=qnorm(0.025, it$K, it$rsCov['K', 'K']), col='red', lty=2)
##
points(skDat$activeCases, col='darkblue', pch=19)
sk$plotMean(add=T, col='blue')
sk$plotBand(col='blue')
#abline(h=sk$K, col='blue')
#abline(h=qnorm(0.975, sk$K, sk$rsCov['K', 'K']), col='blue', lty=3)
#abline(h=qnorm(0.025, sk$K, sk$rsCov['K', 'K']), col='blue', lty=3)
#
points(usDat$activeCases, col='black', pch=19)
us$plotMean(add=T, col='black')
us$plotBand()
#
legend('topleft', legend=c('South Korea', 'Italy', 'USA'), lwd=3, col=c('blue', 'red', 'black'))
dev.off()


