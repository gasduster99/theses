rm(list=ls())

#
library(GA)
library(crch)
library(pracma)
library(EnvStats)
library(numDeriv)
library(rootSolve)
library(latex2exp)
library(RColorBrewer)
#
source('prodClass0.1.1.r')


#
#FUNCTIONS
#

#
dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
	FUpper = exp(lalpha)-M
        Fmsy = uniroot.all(function(ff){ a(gamma, ff, M)-exp(lalpha) }, c(0, FUpper))
        #
        R = exp(lalpha)*P*(1-exp(lbeta)*gamma*P)^(1/gamma) 
        out = R - (M+Fmsy*C)*P
        #
        return( list(out) )
}

#
SRR = function(B, alpha, beta, gamma){
	alpha*B*(1-beta*gamma*B)^(1/gamma)
}

#
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) } #((alpha/(M+ff))^(1/gamma)-1)/beta }

#
getBeta = function(alpha, gamma, M, B0){
	(1-(M/alpha)^gamma)/B0/gamma
}
b = Vectorize(getBeta, c("alpha", "gamma"))

#
getAlpha = function(gamma, ff, M){
	(ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}
a = Vectorize(getAlpha, "gamma")

#
getZeta = function(gamma, ff, M){
	#fgfm = ff*gamma/(ff+M)
	#fgfm/(1 + fgfm - (M/ff+M)^gamma)
	
	#
	(1-((M+ff)/getAlpha(gamma, ff, M))^gamma) / (1-(M/getAlpha(gamma, ff, M))^gamma)
}
z = Vectorize(getZeta, "gamma")

#
getGamma = function(zeta, ff, M){
	uniroot.all(function(gamma){z(gamma, ff, M)-zeta}, c(-500, 100))
	#uniroot.all(function(gamma){ ((ff/(M+ff))*(1-zeta)/(zeta*gamma)+1)^gamma - (M+ff)/M }, c(-1, 100))
	#strongRoot()
}

#
FMsy = function(alpha, gamma, M){
        #
        FUpper = alpha-M
        root = uniroot.all(function(ff){ a(gamma, ff, M)-alpha }, c(0, FUpper))
	
	#root = uniroot.all(function(ff){ 1 - ff/gamma/(M+ff) - (alpha/(M+ff))^(-1/gamma) }, c(0, FUpper))
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        ##
        return(root)
}

#
dttV = function(x, loc, scale, df, V, pV, left=-Inf, right=Inf){
	#out = V*dbinom(n, 1, pV) 
	(1-pV)*dtt(x, loc, scale, df, left, right) #+ pV
}
dttV = Vectorize(dttV, 'x')

#
rttV = function(n, loc, scale, df, V, pV, left=-Inf, right=Inf){
	out = V*rbinom(n, 1, pV) 
	out[out==0] = rtt(sum(out==0), loc, scale, df, left, right)
}
#
pttV = function(x, loc, scale, df, V, pV, left=-Inf, right=Inf){
        #
	out = 0
	if(x>=V){ out=pV + (1-pV)*ptt(x, loc, scale, df, left, right) }
	#
	return(out)
}
pttV = Vectorize(pttV, 'x')

#
#MAIN
#

#
B0 = 10000
M  = 0.2 #05

#
n = 150
xi = 0.5 #3 #0.5 #seq(0, 4, length.out=n+1)
ff = xi*M

#
broke = suppressWarnings(stats::optimize(function(x){z(x, ff, M)}, c(-100, 100)))
minG = broke$minimum
minZ = broke$objective
#
f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy    
fv = Vectorize(f, "x")
peak = stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))
peakG = peak$minimum

#
varZ = 1/optimHess(peakG, function(x){-log(fv(x))})

#
q = 0.9#9

#
l2_dist <- function(f, g) {
  f_diff <- function(x) (f(x) - g(x))^2       
  sqrt(integrate(f = f_diff, lower=minG, upper=qnorm(q, peakG, sqrt(varZ)), subdivisions=10^5)$value)
}

##
#par = c(peakG, sqrt(varZ), 2, minG)
#o = optim(par, 
#	function(x){
#	a = x[1]
#	b = x[2]
#	c = x[3]
#	d = x[4]
#        #(ptt(minG, peakG, x, df)-minZ)^2 +
#        #(ptt(peakG, peakG, x, df)+ptt(minG, peakG, x, df)-z(peakG, ff, M))^2
#        l2_dist(
#                function(y){
#                        ptt(y, a, b, c, left=d)
#                        #ptt(y, a, b, c, left=minG)/(1+ptt(minG, a, b, c))+ptt(minG, a, b, c)
#                },
#                function(y){
#                        z(y, ff, M)
#                }
#        )
#	}
#)
#peakGo = o$par[1]
#sdZ = o$par[2]
#df = o$par[3]
#lef = o$par[4]

#
par = c(peakG, sqrt(varZ), 2)
o = optim(par, 
	function(x){
	#
	a = x[1]
	b = x[2]
	c = x[3]
	#d = x[4]
	#
        l2_dist(#          x, loc,scale,df, V, pV,   left=-Inf, right=Inf
                function(y){#                  d
			pttV(y, a, b, c, minG, minZ, left=minG)
                },
                function(y){
                        z(y, ff, M)
                }
        )
	}
)
peakGo = o$par[1]
sdZo = o$par[2]
dfo = o$par[3]
pVo =  minZ #o1$par[4]

#
par = c(peakG, sqrt(varZ), 2)
o1 = optim(par[-3], 
	function(x){
	#
	a = x[1]
	b = x[2]
	#c = x[3]
	#d = x[4]
	#
        l2_dist(
                function(y){#                  d
			pttV(y, a, b, 100, minG, minZ, left=minG)
                },
                function(y){
                        z(y, ff, M)
                }
        )
	}
)
peakGo1 = o1$par[1]
sdZo1 = o1$par[2]
dfo1 = 100 #o$par[3]
pVo1 =  minZ

#
qq = q #0.99

#
png('zeta.png')
curve(z(x, ff, M), -5, qtt(qq, peakG, sdZo, 2), 
	n=10000, 
	ylim=c(0,1), 
	lwd=3,
	main=TeX("$\\zeta(\\gamma)$"),
	ylab=TeX("$B_{MSY}/B_0$"),
	xlab=TeX("$\\gamma$")
)
points(minG, minZ, pch=19)
legend('bottomleft', legend=TeX('$(\\gamma_{min}, \\zeta_{min})$'), pch=19)
dev.off()


curve(z(x, ff, M), -5, qtt(qq, peakG, sdZo, 2), n=10000, ylim=c(0,1), lwd = 3)
#abline(v=peakG)

#
f = function(x){ pttV(x, peakGo1, sdZo1, dfo1, minG, pVo1, left=minG) }
curve(f(x), minG, qtt(qq, peakG, sdZo, 2), col='green', add=T)

#
g = function(x){ ptt(x, peakG, sqrt(varZ), 2, left=minG) }
curve(g(x), minG, qtt(qq, peakG, sdZo, 2), col='blue', add=T)

#
h = function(x){ pttV(x, peakGo, sdZo, dfo, minG, pVo, left=minG) }
curve(h(x), minG, qtt(qq, peakG, sdZo, 2), col='red', add=T)

#
#FOR THE PRES
#

#
n = 150
xi = 2.5 #3 #0.5 #seq(0, 4, length.out=n+1)
ff = xi*M

#
broke = suppressWarnings(stats::optimize(function(x){z(x, ff, M)}, c(-100, 100)))
minG = broke$minimum
minZ = broke$objective
#
f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy    
fv = Vectorize(f, "x")
peak = stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))
peakG = peak$minimum

#
varZ = 1/optimHess(peakG, function(x){-log(fv(x))})

#
q = 0.9#9

#
l2_dist <- function(f, g) {
  f_diff <- function(x) (f(x) - g(x))^2       
  sqrt(integrate(f = f_diff, lower=minG, upper=qnorm(q, peakG, sqrt(varZ)), subdivisions=10^5)$value)
}

##
#par = c(peakG, sqrt(varZ), 2, minG)
#o = optim(par, 
#	function(x){
#	a = x[1]
#	b = x[2]
#	c = x[3]
#	d = x[4]
#        #(ptt(minG, peakG, x, df)-minZ)^2 +
#        #(ptt(peakG, peakG, x, df)+ptt(minG, peakG, x, df)-z(peakG, ff, M))^2
#        l2_dist(
#                function(y){
#                        ptt(y, a, b, c, left=d)
#                        #ptt(y, a, b, c, left=minG)/(1+ptt(minG, a, b, c))+ptt(minG, a, b, c)
#                },
#                function(y){
#                        z(y, ff, M)
#                }
#        )
#	}
#)
#peakGo = o$par[1]
#sdZ = o$par[2]
#df = o$par[3]
#lef = o$par[4]

#
par = c(peakG, sqrt(varZ), 2)
o = optim(par, 
	function(x){
	#
	a = x[1]
	b = x[2]
	c = x[3]
	#d = x[4]
	#
        l2_dist(#          x, loc,scale,df, V, pV,   left=-Inf, right=Inf
                function(y){#                  d
			pttV(y, a, b, c, minG, minZ, left=minG)
                },
                function(y){
                        z(y, ff, M)
                }
        )
	}
)
peakGo = o$par[1]
sdZo = o$par[2]
dfo = o$par[3]
pVo =  minZ #o1$par[4]



#from=-5
png('zetaTitle.png')
curve(z(x, ff, M), -2, qtt(qq, peakG, sdZo, 2),
        n=10000,
        ylim=c(0,1),
        lwd=3,
        main=TeX("$\\frac{B_{MSY}}{B_0}(\\gamma)$"),
        ylab=TeX("$B_{MSY}/B_0$"),
        xlab=TeX("$\\gamma$")
)
points(minG, minZ, pch=19, cex=1.5)
legend('topright', legend=TeX('$(\\gamma_{min}, \\zeta_{min})$'), pch=19)
dev.off()

#
cols = brewer.pal(9, 'Set1')
png('zetaPlusPlus.png')
curve(z(x, ff, M), -2, qtt(qq, peakG, sdZo, 2),
        n=10000,
        ylim=c(0,1),
        lwd=3,
        main=TeX("$\\frac{B_{MSY}}{B_0}(\\gamma)$"),
        ylab=TeX("$B_{MSY}/B_0$"),
        xlab=TeX("$\\gamma$")
)
points(minG, minZ, pch=19, cex=1.5)
legend('topright', legend=TeX('$(\\gamma_{min}, \\zeta_{min})$'), pch=19)
gs = rttV(150, peakGo, sdZo, dfo, minG, pVo, left=minG)
rug(gs, lwd=1)
rug(z(gs, ff, M), side=2, lwd=1)
curve(dttV(x, peakGo, sdZo, dfo, minG, pVo, left=minG), minG, qtt(qq, peakG, sdZo, 2), add=T, lwd=3, col=cols[2])
segments(minG, 0, minG, minZ, col=cols[2])
points(minG, minZ, pch=20, col=cols[2], cex=1.5)
dev.off()





#
#quantiles
#

#
gs = rttV(10000, peakGo, sdZo, dfo, minG, pVo, left=minG)
zs = z(gs, ff, M)
#
gsBase = rtt(10000, peakG, sqrt(varZ), 2, left=minG)
zsBase = z(gsBase, ff, M)
#
gso1 = rttV(10000, peakGo1, sdZo1, dfo1, minG, pVo1, left=minG)
zso1 = z(gso1, ff, M)

png("qqUnif.png")
qqPlot(zs, 
	distribution="unif", 
	param.list=list(min=minZ, max=1), 
	#add.line=T,
	ylab="Sample Quantiles",
	xlab="Uniform Quantiles",
	main=TeX("Uniform Q-Q Plot for Modeled Samples of $\\zeta$")
)
abline(0, 1, col='red')
dev.off()

##
#dev.new()
#qqPlot(zs, distribution="unif", param.list=list(min=minZ, max=1), add.line=T)
#dev.new()
#qqPlot(zsBase, distribution="unif", param.list=list(min=minZ, max=1), add.line=T)
#dev.new()
#qqPlot(zso1, distribution="unif", param.list=list(min=minZ, max=1), add.line=T)






#dev.new()
#d = function(x){dtt(x, peakG, sqrt(varZ), df, left=minG)}
#curve(d, -5, peakG+qt(0.9, df)*sqrt(varZ))

#
#       minZ = suppressWarnings(stats::optimize(function(x){z(x, ff, M)}, c(-100, 100))$minimum)
#       f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy    
#       fv = Vectorize(f, "x")
#       peakZ = stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))$minimum
#       varZ = 1/optimHess(peakZ, function(x){-log(fv(x))})
#       if(varZ<0){ next }
#       #
#       gs = rtt(1, peakZ, sqrt(varZ)*fudge, df, left=minZ)











##
#out = stats::optimize(function(x){ 
#	#(ptt(minG, peakG, x, df)-minZ)^2 +
#	#(ptt(peakG, peakG, x, df)+ptt(minG, peakG, x, df)-z(peakG, ff, M))^2
#	l2_dist(
#		function(y){
#			ptt(y, peakG, x, df)
#			#ptt(y, peakG, x, df, left=minG)/(1+ptt(minG, peakG, x, df))+ptt(minG, peakG, x, df)
#		},
#		function(y){
#			z(y, ff, M)
#		}
#	)
#	}, c(0, 100))
#sdZ = out$minimum




##number of samples in design
#n = 150 #500 maxed out method
#
##LHS boundaries
#xiMin = 0.25 #3.25 #0.5
#xiMax = 3.75 #3.5 
#zetaMin = 0.15 #0.32 #0.15 
#zetaMax = 0.7  #0.43 #0.7 
##zeta range should not be too small
#
##xi and zeta Bins (bin defined by left edge right egde not used) [n+1 used to give n left edges]
#xiL = seq(xiMin, xiMax, length.out=n+1)
#zetaL = seq(zetaMin, zetaMax, length.out=n+1)
#minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
#binTrk = ceiling(abs(log10(minDiff)))
##
#Llist = rep(T, n) #a list to track which zeta bin still needs a sample
#zetaList = rep(NA, n) #a list the sampled zeta value in each bin
#xiList = rep(NA, n)   #a list the sampled xi value in each bin
#
##details of how to sample the curvy bit for comuting sd
#scope = 2 #how far to sample the inverse function for curvature
#scopeN = 30 #how many samples to take of the inverse function curvature
#scopeRun = seq(-scope, scope, length.out=scopeN)
##start at fudge (hessian sd) slowing increase fudge by fudgeFactor so that we eventually search flat region, but still allow for repetition at calculated sd 
#fudge = 1
#giveUp = 1000
#fudgeFactor = 0.01
##start df at dfStart with thin tails (to explore curvature), decrement by dfFactor to transition to explore flat regiem 
#dfStart = 100
#df = dfStart
#dfFactor = 1 #dfFactor, dfStart, and fudgeFactor are balanced to first decrease df, and then crank on sd
#
##
#xii = 1
#xiShuffle = sample(xiL[-(n+1)])
##
#while( any(Llist) ){
#	#for each xi column (n of them)
#	xi = runif(1, xiShuffle[xii], xiL[which(xiL==xiShuffle[xii])+1])
#	ff = xi*M
#	
#	#
#	minZ = suppressWarnings(stats::optimize(function(x){z(x, ff, M)}, c(-100, 100))$minimum)
#	f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy    
#	fv = Vectorize(f, "x")
#	peakZ = stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))$minimum
#	varZ = 1/optimHess(peakZ, function(x){-log(fv(x))})
#	if(varZ<0){ next }
#	#
#	gs = rtt(1, peakZ, sqrt(varZ)*fudge, df, left=minZ)
#	
#	#propose sample of zeta
#	zs = z(gs, ff, M)
#	#reject: catch the case where zs is propsed outside of z range
#	if(zs<zetaMin | zs>zetaMax){ next } 
#	
#	#NOTE: compute auxilliary quantities
#	#NOTE: make sure that reverse a,b,g->xi,zeta and xi,zeta->a,b,g calcs make sense; else reject	
#	as = a(gs, ff, M)
#	bs = b(as, gs, M, B0)
#	
#	#
#	fStar = FMsy(as, gs, M)
#	if( length(fStar)==0 ){ 
#		#print('Skip Fmsy')
#		next
#	}
#	#BStar = PBar(as, bs, gs, ff, M)
#	#BZero = PBar(as, bs, gs, 0, M)	
#	##print(c(as, bs, gs))
#	#print(c(xi, FMsy(as, gs, M)/M))
#	#print(c(ff, FMsy(as, gs, M)))
#	#print(c(zs, BStar/BZero))
#	#print('')
#	##if( reverse a,b,g->xi,zeta and xi,zeta->a,b,g fail ){ next } #reject
#
#	#find left edge of bin sampled
#	bin = max(which(zetaL<=zs))
#	#if sampled bin not yet sampled, then save
#	if( Llist[bin] ){
#		#record values
#		xiList[bin] = xi 
#		zetaList[bin] = zs
#		#
#		#NOTE: generate a faux data object and catelogue design elements in a model object in the design directory
#		#make Schunte production model
#               	datGen = prodModel$new(
#               	        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
#               	        time=1:30, catch=FtFmsy, M=M,                   #constants
#               	        alpha=as, beta=bs, gamma=gs,                   #parameters
#               	        lalpha=log(as), lbeta=log(bs),                 #reparameterize
#               	        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
#               	        xi=xi, zeta=zs                                #other incidentals to carr
#               	)
#               	datGen$iterate(odeMethod)
#		datGen$plotMean()
#	
#		#NOTE: catelogue by bin lables (left edges of xiL, zetaL)s  xiShuffle[xii]
#		datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiShuffle[xii], binTrk), round(zetaL[bin], binTrk)))
#		#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiL[bin], binTrk), round(zetaL[bin], binTrk)))
#
#		#knobs for making the design
#		Llist[bin] = F 	#turn off bin
#		xii = xii+1 	#increment xi counter
#		df = 100    	#reset df to thin tails
#		fudge = 1   	#reset fudge to thin tails
#	}else{ # catch case where hyperbola gap requires heavy tails
#		df = max(1, df-1)  #decrement df to broadend search to get the center bald spot
#		fudge = fudge+0.01 #increase fudge factor to get a wider sd
#		if( fudge>giveUp ){ break } #give up (it mostly just misses one of two bins so its usually fine)
#	}
#}

#
#JUNK
#

##
#place = './modsSchnuteHHardExpT45N150M0.05Wide/'
#odeMethod = "lsode"
#
##safely create place
#if( dir.exists(place) ){
#	#
#	isGo = readline(sprintf("Overwrite %s Design?\nCtrl-C to Escape, Enter to Overwrite, or A to Append: ", place))	
#	#
#	if( toupper(isGo)=="A" ){
#		#
#		dir.create(place, showWarnings=FALSE)
#	}else{
#		#
#		unlink(place, recursive=TRUE)
#		dir.create(place, showWarnings=FALSE)
#	}
#} else{ dir.create(place, showWarnings=FALSE) }
#
##
#TT = 30
#FtFmsy = rep(1, TT) #make faux catch



##
##Extreme Design Elements
##
#
##c(xiMax, zetaMax)
#xi   = xiMax
#zeta = zetaMax
##curve(z(x, xi*M, M), 0, 10)
#gs = uniroot.all(function(gamma){z(gamma, xi*M, M)-zeta}, c(0, 10))
#as = a(gs, xi*M, M)
#bs = b(as, gs, M, B0)
##make Schunte production model
#datGen = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
#        time=1:30, catch=FtFmsy, M=M,                  #constants
#        alpha=as, beta=bs, gamma=gs,                   #parameters
#        lalpha=log(as), lbeta=log(bs),                 #reparameterize
#        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
#        xi=xi, zeta=zeta                               #other incidentals to carr
#)
#datGen$iterate(odeMethod)
#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xi, binTrk), round(zeta, binTrk)))
#
##c(xiMax, zetaMin)
#xi   = xiMax
#zeta = zetaMin
##curve(z(x, xi*M, M), -1, 10)
#gs = uniroot.all(function(gamma){z(gamma, xi*M, M)-zeta}, c(-1, 10))
#as = a(gs, xi*M, M)
#bs = b(as, gs, M, B0)
##make Schunte production model
#datGen = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
#        time=1:30, catch=FtFmsy, M=M,                  #constants
#        alpha=as, beta=bs, gamma=gs,                   #parameters
#        lalpha=log(as), lbeta=log(bs),                 #reparameterize
#        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
#        xi=xi, zeta=zeta                               #other incidentals to carr
#)
#datGen$iterate(odeMethod)
#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xi, binTrk), round(zeta, binTrk)))
#
##c(xiMin, zetaMax)
#xi   = xiMin
#zeta = zetaMax
##curve(z(x, xi*M, M), -1, 10)
#gs = uniroot.all(function(gamma){z(gamma, xi*M, M)-zeta}, c(-1, 10))
#as = a(gs, xi*M, M)
#bs = b(as, gs, M, B0)
##make Schunte production model
#datGen = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
#        time=1:30, catch=FtFmsy, M=M,                  #constants
#        alpha=as, beta=bs, gamma=gs,                   #parameters
#        lalpha=log(as), lbeta=log(bs),                 #reparameterize
#        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
#        xi=xi, zeta=zeta                               #other incidentals to carr
#)
#datGen$iterate(odeMethod)
#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xi, binTrk), round(zeta, binTrk)))
#
##c(xiMin, z(xiMin))
#xi = xiMin
#zetaOpt = suppressWarnings(stats::optimize(function(x){z(x, xi*M, M)}, c(-100, 100)))
#zeta = zetaOpt$objective #+ 2*.Machine$double.eps^0.25
##curve(z(x, xi*M, M), -1, 10)
#gs = zetaOpt$minimum 	#+ 2*.Machine$double.eps^0.25
#as = 0.969 #1 			#a(gs, xi*M, M)
#bs = b(as, gs, M, B0)
##make Schunte production model
#datGen = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
#        time=1:30, catch=FtFmsy, M=M,                  #constants
#        alpha=as, beta=bs, gamma=gs,                   #parameters
#        lalpha=log(as), lbeta=log(bs),                 #reparameterize
#        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
#        xi=FMsy(as, gs, M)/M, zeta=getZeta(gs, FMsy(1, gs, M), M)                               #other incidentals to carr
#)
#datGen$iterate(odeMethod)
#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xi, binTrk), round(zeta, binTrk)))
#
###c(xiMax, z(xiMax))
##xi = xiMax
##zetaOpt = suppressWarnings(stats::optimize(function(x){z(x, xi*M, M)}, c(-100, 100)))
##zeta = zetaOpt$objective #+ 2*.Machine$double.eps^0.25
###curve(z(x, xi*M, M), -1, 10)
##gs = zetaOpt$minimum 	#+ 2*.Machine$double.eps^0.25
##as = a(gs, xi*M, M)
##bs = b(as, gs, M, B0)
###make Schunte production model
##datGen = prodModel$new(
##        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
##        time=1:30, catch=FtFmsy, M=M,                  #constants
##        alpha=as, beta=bs, gamma=gs,                   #parameters
##        lalpha=log(as), lbeta=log(bs),                 #reparameterize
##        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
##        xi=FMsy(as, gs, M)/M, zeta=getZeta(gs, FMsy(1, gs, M), M)                               #other incidentals to carr
##)
##datGen$iterate(odeMethod)
##datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xi, binTrk), round(zeta, binTrk)))








##
#png('design.png')
#plot(xiList, zetaList, pch=20)
##abline(v=xiL, lty=3)
##abline(h=zetaL, lty=3)
#curve(1/(x+2), 0, 4, add=T)
#dev.off()

##
#f = function(x){ ((ff/(M+ff))*(1-zeta)/(zeta*x)+1)^x - (M+ff)/M }
##uniroot(f, c(-1, -zeta-0.1))

##
#ff=0.2
#gamma = 1
#
#
#zz = getZeta(gamma, ff, M)
#alpha = getAlpha(gamma, ff, M)
#beta = getBeta(alpha, gamma, M, B0)
###print( PBar(alpha, beta, gamma, 0, M) )
##
###
#plot(SRR(0:B0, alpha, beta, gamma)-0:B0*M, type='l')
##curve(SRR(x, alpha, beta, gamma), from=0, to=B0)
#lines(0:B0, 0:B0*(ff))#+M))
#abline(v=zeta*B0)

## 
	#pole = getPole(ff, M)
	#
	##tune scale on curvature
	#f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy	
	#fv = Vectorize(f, "x")
	#hesses = sapply(scopeRun+pole, function(i){numDeriv::hessian(function(x){log(f(x))}, i)})	
	#sd = median(sqrt(-1/hesses), na.rm=T)
	##may (maybe not) be beneficial to consider adding a tuning step for the flatter sections of getZeta(z)??
	#
	#cub = theWolf(ff, M)
	##tune df on student t tails and thus emphasize or de-emphasize the central peak.
	#smaller df has wider tails and lower spread (since it samples flat regiem of getZeta), 
	#higher df has thinner tails and more spread (since is samples the spiky bit of getZeta). 
	#gs = cub$pole + rt(1, df)*cub$sd*fudge #rt(1, 100)*sd*fudge  
	#wolfSample(ff, M, df, fudge)
	
##
#getPole = function(ff, M){
#	##
#	#if(ff<=M){
#	#	out = max(uniroot.all( 
#	#		function(gamma){
#	#			fgfm = ff*gamma/(ff+M) 
#	#			(1 + fgfm - (M/ff+M)^gamma)
#	#		},  c(eps(),100) 
#	#	))
#	#} else if(ff>M){
#		f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy     
#        	fv = Vectorize(f, "x")
#		out = optimize(function(x){-fv(x)}, c(-10, 100))$minimum
#	#}
#	
#	#
#	return(out)
#}

##
#wolfSample = function(ff, M, df, fudge){
#	#tune scale on curvature
#        f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy     
#        fv = Vectorize(f, "x")
#	pole = optimize(function(x){-fv(x)}, c(-100, 100))$minimum
#	
#	#
#        if(ff<=M){
#               	##
#		#pole = max(uniroot.all( 
#                #       function(gamma){
#                #               fgfm = ff*gamma/(ff+M) 
#                #               (1 + fgfm - (M/ff+M)^gamma)
#                #       },  c(-100,100) 
#		#))
#		
#		#
#		#pole = optimize(function(x){-fv(x)}, c(-100, 100))$minimum
#		hesses = sapply(scopeRun+pole, function(i){numDeriv::hessian(function(x){log(f(x))}, i)})
#       		sd = median(sqrt(1/hesses), na.rm=T)
#		#
#		gs = pole + -abs(rt(1, df)*sd*fudge)
#        } else if(ff>M){
#                #pole = optimize(function(x){-fv(x)}, c(-100, 100))$minimum
#        	
#		#
#		hesses = sapply(scopeRun+pole, function(i){numDeriv::hessian(function(x){log(f(x))}, i)})
#        	sd = median(sqrt(-1/hesses), na.rm=T)
#		#
#		gs = pole + rt(1, df)*sd*fudge
#	}
#	#
#	return(gs)
#}
##
#xi = runif(1, xiShuffle[xii], xiL[which(xiL==xiShuffle[xii])+1])
#ff = xi*M
##
#minZ = suppressWarnings(optimize(function(x){z(x, ff, M)}, c(-20, 100))$minimum)
#f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy    
#fv = Vectorize(f, "x")
#peakZ = optimize(function(x){-log(fv(x))}, c(minZ, 100))$minimum
#varZ = 1/optimHess(peakZ, function(x){-log(fv(x))})
##
#curve(fv(x), minZ, 20, ylim=c(0, 1), n=1000)
#curve(dtt(x, peakZ, sqrt(varZ), 100, left=minZ), minZ, 20, add=T, col='red')


