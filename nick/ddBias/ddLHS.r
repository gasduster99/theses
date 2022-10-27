rm(list=ls())

#
library(crch)
library(Ryacas)
library(pracma)
library(stringr)
library(rootSolve)

#
source('ddClass0.1.1.r')

#
#FUNCTIONS
#

#prepmath
##
#R = ysym("alpha*B*(1-gamma*beta*B)^(1/gamma)")
##
#Nbar = ysym("R - (M+F)*N")
#Nbar = with_value(Nbar, "R", R)
#Nbar = ysym("(alpha*B*(1-gamma*beta*B)^(1/gamma))/(M+F)") #solve(Nbar, 'N')[1]
##
##Bbar = ysym("w*R + k*(W*N-B) - (M+F)*B")
##Bbar = with_value(Bbar, "R", R)
##Bbar = with_value(Bbar, "N", Nbar)

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%28W*
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
Bbar_r = as_r( with_value(Bbar, "F", ysym("FF")) )
BBar = function(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#
FBbar = ysym("F*B")
FBbar = with_value(FBbar, "B", Bbar)
FBbar = with_value(FBbar, "F", ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
dFBdF_r = as_r(dFBdF)
FMsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10))$root
}
#beta does not matter for either of these
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma, dFBdFexpr=dFBdF_r){
	#
	#capture.output(
	out <- tryCatch({
        	uniroot(function(alpha){ eval(dFBdFexpr) }, c(eps(), 100))$root
	}, error=function(err){
		out = NA
	})#, file="/dev/null")
	#
	return(out)
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta, dFBdFexpr=dFBdF_r){
        #FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(dFBdFexpr) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
        f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
        uniroot(f, c(0, 10), tol=eps())$root
}

#gamma|alpha, zeta
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
        f = function(g){
                BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - zeta
        }
        uniroot(f, c(-10, 10))$root
}

#
getZeta = function(FF, M, k, w, W, alpha, beta, gamma){
	BBar(FF, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
}

#
invert = function(zeta, xi, B0, M, k, w, W, alphaStart=1, betaStart=1, gammaStart=1){
	#
	alpha = alphaStart
	beta  = betaStart
	gamma = gammaStart
	#
	xTol = 10^-3
	zTol = xTol
	bTol = 10
	xiHat = xi+xTol*10
	zetaHat = zeta+zTol*10
	#
	while( abs(xiHat-xi)>=xTol | abs(zetaHat-zeta)>=zTol ){
	        #
	        gamma = getGammaZeta(zeta, xi*M, M, k, w, W, alpha, beta)
	        alpha = getAlphaFmsy(xi*M, M, k, w, W, beta, gamma)
	        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
	        #
	        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
	        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
	        zetaHat = BBar(xiHat*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
	        #
	        #xis = c(xis, xiHat)
	        #zetas = c(zetas, zetaHat)
	        #i = i+1
	}
	#
	return( data.frame(alpha=alpha, beta=beta, gamma=gamma) )
}


#
der = function(t, Y, lalpha, lbeta, gamma, a0, WW, kappa, catch, B0){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        N = Y[1]
        B = Y[2]
        #
        if( (t-a0)<1){
                Blag = B0
        }else{
                Blag = lagvalue(t-a0)[2]
        }
        #
        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        alpha = exp(lalpha)
        beta = exp(lbeta)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
        #
	ww = WW*(1-exp(-kappa*a0))
	FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
        #
        out = c(N=NA, B=NA)
        out[1] = R - (M+FF)*N
        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B 
        #
        return( list(out) )
}

#
#MAIN
#

#
a0 = 2
M  = 0.2
kappa = 0.2
WW = 1
ww = WW*(1-exp(-kappa*a0))

#
B0 = 10000

#
xi = 1
zeta = 0.4

#
inv = invert(zeta, xi, B0, M, kappa, ww, WW)
alpha = inv$alpha
beta  = inv$beta
gamma = inv$gamma

#
TT = 45
FtFmsy = rep(1, TT) #make faux catch

#
Fmsy  = FMsy(M, kappa, ww, WW, alpha, beta, gamma)
BStar = BBar(Fmsy, M, kappa, ww, WW, alpha, beta, gamma)
BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)

#
#LHS
#

#BBar goes negative if ff or xi gets too large
#ffTop = uniroot(function(x)BBar(x, M, kappa, ww, WW, alpha, beta, gamma), c(0.4, 0.5))$root
xiLim = c(0.5, 3.75) #c(0.5, ffTop/M)  #c(0.5, 2.4) #c(0.5, 3.75)
zetaLim = c(0.15, 0.7)
n = 100

#LHS boundaries
xiMin = xiLim[1]        #0.25 #3.25 #0.5
xiMax = xiLim[2]        #3.75 #3.5 
zetaMin = zetaLim[1]    #0.15 #0.32 #0.15 
zetaMax = zetaLim[2]    #0.7  #0.43 #0.7 
#zeta range should not be too small

#xi and zeta Bins (bin defined by left edge right egde not used) [n+1 used to give n left 
xiL = seq(xiMin, xiMax, length.out=n+1)
zetaL = seq(zetaMin, zetaMax, length.out=n+1)
minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
binTrk = ceiling(abs(log10(minDiff)))
#
Llist = rep(T, n) #a list to track which zeta bin still needs a sample
zetaList = rep(NA, n) #a list the sampled zeta value in each bin
xiList = rep(NA, n)   #a list the sampled xi value in each bin

#details of how to sample the curvy bit for comuting sd
scope = 2 #how far to sample the inverse function for curvature
scopeN = 30 #how many samples to take of the inverse function curvature
scopeRun = seq(-scope, scope, length.out=scopeN)
#start at fudge (hessian sd) slowing increase fudge by fudgeFactor so that we eventually s
fudge = 1
giveUp = 100 #10
fudgeFactor = 0.1 #0.01
#start df at dfStart with thin tails (to explore curvature), decrement by dfFactor to tran
dfStart = 100
df = dfStart
dfFactor = 1 #dfFactor, dfStart, and fudgeFactor are balanced to first decrease df, and th

#
xii = 1
xiShuffle = sample(xiL[-(n+1)])
#
while( any(Llist) ){
	#for each xi column (n of them)
	xi = runif(1, xiShuffle[xii], xiL[which(xiL==xiShuffle[xii])+1])
	ff = xi*M
	
	#should return NA if something fails so that minZ works as planned
	z = function(g){ 
		a = getAlphaFmsy(ff, M, kappa, ww, WW, beta, g)
		getZeta(ff, M, kappa, ww, WW, a, beta, g) 
	}
	z = Vectorize(z, "g")
	minZ = suppressWarnings(stats::optimize(z, c(-100, 100))$minimum)
	f = function(x){grad(z, x)}#pdf analogy    
	fv = Vectorize(f, "x")
	peakZ = suppressWarnings(stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))$minimum)
	varZ = tryCatch({ 1/optimHess(peakZ, function(x){-log(fv(x))})
	}, error=function(err){ varZ=-1 })
	if(varZ<0){ next }
        #
	gs = rtt(1, peakZ, sqrt(varZ)*fudge, df, left=minZ)

	#propose sample of zeta
	zs = z(gs)
	#reject: catch the case where zs is propsed outside of z range
	if(zs<zetaMin | zs>zetaMax){ next }#NOTE: if slow, I could probably use the invert function to include bounds in truncation above	
	
	#
	as = getAlphaFmsy(ff, M, kappa, ww, WW, beta, gs)
	bs = getBeta(B0, M, kappa, ww, WW, as, gs)

	#
	fStar = FMsy(M, kappa, ww, WW, as, bs, gs)
	#if( length(fStar)==0 ){
	#	#print('Skip Fmsy')
	#	next
	#}
	
	#find left edge of bin sampled
	bin = max(which(zetaL<=zs))
	#if sampled bin not yet sampled, then save
	if( Llist[bin] ){
	        #record values
	        xiList[bin] = xi
	        zetaList[bin] = zs
		
		#
		dat = ddModel$new( derivs=der,
			N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){	#virgin numbers
			        #
			        alpha = exp(lalpha)
			        beta  = exp(lbeta)
				ww = WW*(1-exp(-kappa*a0))
			        #
				BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma) 
			        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
				
			},
			B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){	#virgin biomass
			        #
			        alpha = exp(lalpha)
			        beta = exp(lbeta)
				ww = WW*(1-exp(-kappa*a0))
			        #
				BBar(0, M, kappa, ww, WW, alpha, beta, gamma) 	
			},
			time=1:TT, catch=FtFmsy, a0=a0, M=M, WW=WW, kappa=kappa,#constants
			lalpha=log(as), lbeta=log(bs), gamma=gs, 		#parameters
			lq=log(0.00049), lsdo=log(0.01160256)           	#nuisance parameters
		)
		dat$iterate()#"radau")
		dat$plotQuan()
		
		#knobs for making the design
		Llist[bin] = F  #turn off bin
		xii = xii+1     #increment xi counter
		df = 100        #reset df to thin tails
		fudge = 1       #reset fudge to thin tails
	}else{ # catch case where hyperbola gap requires heavy tails
		df = max(1, df-1)  #decrement df to broadend search to get the center bald spot
		fudge = fudge+fudgeFactor #increase fudge factor to get a wider sd
		if( fudge>giveUp ){ break } #give up (it mostly just misses one of two bins so its usually fine)
	}
}


#N0 = ( alpha/(beta*(1+gamma)) * (1-gamma/(1+gamma))^(1/gamma) )/M
#B0 = ( WW*(1-exp(-kappa*a0))*(alpha/(beta*(1+gamma))*(1-gamma/(1+gamma))^(1/gamma)) +
#	kappa*WW*(alpha/(beta*(1+gamma))*(1-gamma/(1+gamma))^(1/gamma))/M
#	)/(kappa+M)

##
#sam = tgp::lhs(200, rbind(xiLim, zetaLim))
#sam = sam[order(sam[,2], decreasing=T),]
#abg = c()
#for(i in 1:nrow(sam)){
#	#
#	inv = invert(sam[i,2], sam[i,1], B0, M, kappa, ww, WW, alphaStart=alpha, betaStart=beta, gammaStart=gamma)
#	alpha = inv$alpha
#	beta  = inv$beta
#	gamma = inv$gamma
#	abg = rbind(abg, c(alpha=alpha, beta=beta, gamma=gamma))
#}



