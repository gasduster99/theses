rm(list=ls())

#
library(plgp)
library(crch)
library(gMOIP)
library(Ryacas)
library(pracma)
library(stringr)
library(rootSolve)

#
source('ddClass0.1.1.r')

#
#FUNCTIONS
#

#EQUILIBRIUM EQUATIONS

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
#
NBar = function(FF, M, k, w, W, alpha, beta, gamma){
        #
        BB = BBar(FF, M, kappa, ww, WW, alpha, beta, gamma)
        (alpha*BB*( 1-beta*gamma*BB )^(1/gamma))/(M+FF)
}
NCar = function(x, FF, M, k, w, W, alpha, beta, gamma){
        #
        #BB = BBar(FF, M, kappa, ww, WW, alpha, beta, gamma)
        (alpha*x*( 1-beta*gamma*x )^(1/gamma))/(M+FF)
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
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#
getZetaBH = function(x, M, k, W, aS, a0){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
        #
        gamma = -1
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaBH = Vectorize(getZetaBH, "x")

#DIFFERENTIAL EQUATION 

#
der = function(t, Y, lalpha, lbeta, gamma, aS, a0, WW, kappa, catch, B0){
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
        if( (t-aS)<1){
                Blag = B0
        }else{
                Blag = lagvalue(t-aS)[2]
        }
        #
        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        alpha = exp(lalpha)
        beta = exp(lbeta)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
        #
	ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
	FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
        #
        out = c(N=NA, B=NA)
        out[1] = R - (M+FF)*N
        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B 
        #
        return( list(out, ww*R, kappa*(WW*N-B)) )
}

#
SRR = function(B, mod){
	#
	alpha = exp(mod$lalpha)
        beta = exp(mod$lbeta)
	gamma = mod$gamma
	#
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
	#
	return(R)
}
SRR = Vectorize(SRR, 'B')

#DESIGN STUFF

#
lhsForPlot = function(xiLim, zetaLim, n, batch, save=F){
        #somehow define place when save!=F
        #deal with append behaviour	
	
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
	mod = vector(mode='list', length=n)
	Llist = rep(T, n) #a list to track which zeta bin still needs a sample
	zetaList = rep(NA, n) #a list of the sampled zeta value in each bin
	xiList = rep(NA, n)   #a list of the sampled xi value in each bin	

	#details of how to sample the curvy bit for comuting sd
	scope = 2 #how far to sample the inverse function for curvature
	scopeN = 30 #how many samples to take of the inverse function curvature
	scopeRun = seq(-scope, scope, length.out=scopeN)
	#start at fudge (hessian sd) slowing increase fudge by fudgeFactor so that we eventually s
	fudge = 1
	giveUp = 100 		#100 #10
	fudgeFactor = 0.1 	#0.1 #0.01
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
		xiLWho = which(xiL==xiShuffle[xii])
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
		if(varZ<0){ 
			#
			fudge = fudge+fudgeFactor #increase fudge factor to make progress toward giving up
                        if( fudge>giveUp ){ break }
			#
			next 
		}
	        #
		gs = rtt(1, peakZ, sqrt(varZ)*fudge, df, left=minZ)
	
		#propose sample of zeta
		zs = z(gs)
		#reject: catch the case where zs is propsed outside of z range
		if(zs<zetaMin | zs>zetaMax){
			# 
			fudge = fudge+fudgeFactor #increase fudge factor to make progress toward giving up
                        if( fudge>giveUp ){ break }
			#
			next 
		}#NOTE: if slow, I could probably use the invert function to include bounds in truncation above	
		
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
					ww = vbGrow(aS, k, W, a0) #WW*(1-exp(-kappa*a0))
				        #
					BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma) 
				        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
					
				},
				B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){	#virgin biomass
				        #
				        alpha = exp(lalpha)
				        beta = exp(lbeta)
					ww = vbGrow(aS, k, W, a0) #WW*(1-exp(-kappa*a0))
				        #
					BBar(0, M, kappa, ww, WW, alpha, beta, gamma) 	
				},
				time=1:TT, catch=FtFmsy, aS=aS, a0=a0, M=M, WW=WW, kappa=kappa,#constants
				lalpha=log(as), lbeta=log(bs), gamma=gs, 		#parameters
				lq=log(0.00049), lsdo=log(0.01160256),	        	#nuisance parameters
				xi=xi, zeta=zs
			)
			dat$iterate(odeMethod)
			mod[[bin]] = dat
			#dat$plotQuan(main=sprintf("%s, %s", xi, zs))
			#writeLines(sprintf("\n\n%s, %s", xi, zs))
			#dat$printSelf()
			
			#       
                        if( save!=F ){
                                if( save==T ){
                                        #create place
                                        place = date()
					save = place
					dir.create(place, showWarnings=FALSE)
                                }else{ place=save }
                                #NOTE: catelogue by bin lables (left edges of xiL, zetaL)s  xiShuffle[xii]
                                #datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiL[bin], binTrk), round(zetaL[bin], binTrk)))
                                datName = sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiShuffle[xii], binTrk), round(zetaL[bin], binTrk))
                                dat$save(datName)
				#
				apPlace = sprintf("%s/appends.csv", place)
				if( !file.exists(apPlace) ){ 
					write.table(t(c('batch', 'modName')), apPlace, sep=',', row.names=F, col.names=F, append=F)
				}
                                write.table(t(c(batch, datName)), apPlace, sep=',', row.names=F, col.names=F, append=T)
                        }
			
			#knobs for making the design
			Llist[bin] = F  #turn off bin
			xii = xii+1     #increment xi counter
			df = 100        #reset df to thin tails
			fudge = 1       #reset fudge to thin tails
		}else{ # catch case where hyperbola gap requires heavy tails
			df = max(1, df-1)  #decrement df to broadend search to get the center bald spot
			fudge = fudge+fudgeFactor #increase fudge factor to get a wider sd
			if( fudge>giveUp ){ break } #give up (it mostly just misses one or two bins so its usually fine)
			#writeLines(sprintf('\t%s', fudge))
		}
	}
	#
	out = list(
		ll = cbind(xiList[!is.na(xiList)], zetaList[!is.na(zetaList)]),
		mod = mod
	)
	#
	return( out )
}

#
move01tolim = function(coord01, lim){
        coord01*(lim[2]-lim[1])+lim[1]
}

#
movelimto01 = function(coordlim, lim){
        (coordlim-lim[1])/(lim[2]-lim[1])
}

#
maximinAddHullScale = function(n, xlim, ylim, TT=1e5, XStart=NULL) {
        #n : number of points to draw
        #m : dimension of space to draw from
        #
        #value :  only returns new additions

        #
        ZStart = cbind( movelimto01(XStart[,1], xlim), movelimto01(XStart[,2], ylim) )
        vchull = ZStart[chull(ZStart),]
        #plot(ZStart)
        #polygon(vchull)

        #
        Z = matrix(NA, nrow=n, ncol=2)
        Z[,1] = rep(vchull[1,1], n) ## start everything in a place that is guaranteed to be inside the hull
        Z[,2] = rep(vchull[1,2], n)
        #
        d = distance(Z)
        d = d[upper.tri(d)]
        #
        md = min(d)
        #handle initiated design
        if(!is.null(ZStart)){
                md2 = min(distance(Z, ZStart))
                if(md2 < md){ md=md2 }
        }
        #
        propZ = cbind(runif(TT), runif(TT))
        isIn = inHull(propZ, vchull)>=0
        inZ = propZ[isIn,]
        #
        tList = 1:sum(isIn)
        for(t in tList){
                #propse a row
                row = sample(1:n, 1)
                zold = Z[row,]
                #assume we accept
                Z[row,] = inZ[t,] ## random new row 
                d = distance(Z) #new distance
                d = d[upper.tri(d)]
                mdprime = min(d)
                #handle initial design (make min dist smaller if neccesary with DStart included)
                if(!is.null(ZStart)){
                        mdprime2 = min(distance(Z, ZStart))
                        if(mdprime2 < mdprime) mdprime=mdprime2
                }
                #
                if( mdprime > md ){ md=mdprime ## accept 
                } else{ Z[row,]=zold } ## reject

        }
        #points(Z, pch=19)
        #addCircle(Z[,1], Z[,2], sqrt(md))

        #
        X = cbind( move01tolim(Z[,1], xlim), move01tolim(Z[,2], ylim) )
        out = list(X=X, Z=Z, md=sqrt(md))

        #
        return( out )
}

#
addCircle = function(centerx, centery, radius, length=200){
        # prepare "circle data"
        theta = seq(0, 2 * pi, length=length) # angles for drawing points around the circle

        # draw the circle
        lines(x = radius * cos(theta) + centerx, y = radius * sin(theta) + centery)
}


#
#HEAD
#

#DATA STUFF

#
#mod = "FlatT45N300A0-1AS10K0.1SDORAND"
#mod = "FlatT45N300A0-1AS10K0.1SDO"
#mod = "FlatT30N150A15K0.1" #"ExpT45N150A15K0.1" #"ExpT45N150A15" #"ExpT45N150K1" #"ExpT45N150A15" # "ExpT45N150Wide" #"ExpT45N150A15K0.1" #"ExpT45N150K1" #
#mod = "ExpT45N150A-1AS15K0.1"
#mod = "ExpT45N300AS0.1K10N56" #"ExpT45N150A-1AS2"
#mod = "ExpT45N300AS10K10" #"ExpT45N300AS1K1N28"
#mod = "ExpT45N300AS0.1K10"
#mod = "ExpT45N300AS10K0.1"
mod = "ExpT45N300A0-1AS10K0.1"
place = sprintf("./modsDD%s/", mod)

#
datFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("datGen_*.rda")))
fitFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("fit_*.rda")))

#
rn = c()
xiList = c()
zetaList = c()
for(fitF in fitFiles){
        #
        fit = readRDS(fitF)
        #check rsCOV
        if(length(fit$rsCov)==0){ next }
        #
	rn = c(rn, gsub(sprintf('%sfit_', place), '', fitF))
        xiList = c(xiList, fit$xi)
        zetaList = c(zetaList, fit$zeta)	
}
l = as.data.frame(cbind(xiList, zetaList))
rownames(l) = rn 

#DD MODEL STUFF

#
odeMethod = "lsode"

#
i = 1
C = NULL
while( is.null(C) ){
        #
        f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda"))[i] )
        fOne = readRDS(f)
        C = fOne$rsCov
        #
        i = i+1
}

#
aS = fOne$aS
a0 = fOne$a0
M  = fOne$M
kappa = fOne$kappa
WW = fOne$WW
ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
#
B0 = 10000

##
#aS = 
#a0 = 15     #15  #7.5 #15  #1
#M  = 0.2
#kappa = 0.1 #0.1 #0.2 #0.2 #10
#WW = 1
#ww = vbGrow(aS, k, W, a0) #WW*(1-exp(-kappa*a0))
##
#B0 = 10000

#just initiating to give alpha, beta, gamma some reasonable values
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

#LAYOUT STUFF

#NOTE:high zeta and high xi lets the initial drop get below Bmsy
N = 3 #4#20

#
xiLim = c(0.5, 3.5)   #c(0.25, 3.75)
zetaLim = c(0.2, 0.65) #c(0.15, 0.7) #0.6) #

#LHS boundaries
xiMin = xiLim[1]        #0.25 #3.25 #0.5
xiMax = xiLim[2]        #3.75 #3.5 
zetaMin = zetaLim[1]    #0.15 #0.32 #0.15 
zetaMax = zetaLim[2]    #0.7  #0.43 #0.7 
#zeta range should not be too small

#xi and zeta Bins (bin defined by left edge right egde not used) [n+1 used to give n left 
xiL = seq(xiMin, xiMax, length.out=N+1)
zetaL = seq(zetaMin, zetaMax, length.out=N+1)

#if sampled bin not yet sampled, then save
l = l[	l$xiList>=xiLim[1] & 
	l$xiList<=xiLim[2] &
	l$zetaList>=zetaLim[1] &
	l$zetaList<=zetaLim[2]
,]
Lay = matrix(0, N, N)
for( i in 1:nrow(l) ){
	#
	xi = l[i,1]
	zs = l[i,2]
        #find left edge of bin sampled
	xBin = max(which(xiL<=xi))
	zBin = max(which(zetaL<=zs))
	#
	Lay[xBin, zBin] = i
}
#
tops = unique(sort(Lay))

#
lay = Lay
j = 1
for(i in 1:nrow(l)){
	#
	if( !i%in%tops ){next}
	#
	lay[lay==i] = j
	j = j+1
}
lay = lay[seq(N,1),] #matrix(lay[seq(N^2,1)], N, N) #
#anti-transpose
P = diag(N)[seq(N,1),]
lay = P%*%t(lay)%*%P

#PLOT

#
png(sprintf('indexGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit = readRDS(fWho)
	#print(fWho)
	#print( FMsy(M, fit$kappa, vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0), fit$WW, fit$alpha, fit$beta, fit$gamma) )
	fit$printSelf()
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat = readRDS(dWho)
	#
	par(mar=c(1,1,1,0))
	dat$plotMean() #main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,dat$B0) )
	fit$plotMean( add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
	fit$plotBand( col='red' )
}
dev.off()

#
il = c()
ll = c()
png(sprintf('indexGridKA%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
        #
        if( !i%in%tops ){next}
	#
	il = c(il, i)
	ll = c(ll, rownames(l)[i])
	#
        fWho = sprintf('%sfit_%s', place, rownames(l)[i])
        fit = readRDS(fWho)
        #
	fBHKAWho = sprintf('%sfitBHKA_%s', place, rownames(l)[i])
        fitBHKA = readRDS(fBHKAWho)
	
	#print(fWho)
        #print( FMsy(M, fit$kappa, vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0), fit$WW, fit$alpha, fit$beta, fit$gamma) )
        #fit$printSelf()
        
	#
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
        dat = readRDS(dWho)
	cpue = rlnorm(length(dat$B), dat$lq+log(dat$B), exp(dat$lsdo)*10)
        #
        par(mar=c(1,1,1,0))
        plot(cpue, cex=2, pch=19)
	#plot(fitBHKA$cpue, cex=2, pch=19)
	#dat$plotMean() #main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,dat$B0) )
        #fit$plotMean( add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$
        #fit$plotBand( col='red' )
	fitBHKA$plotMean( add=T, col='blue')
	fitBHKA$plotBand( col='blue' )
	
	#
	f3KAWho = sprintf('%sfit3KA_%s', place, rownames(l)[i])
        if( file.exists(f3KAWho) ){ 
		fit3KA = readRDS(f3KAWho)
		fit3KA$plotMean( add=T, col='green')
        	fit3KA$plotBand( col='green' )
	}else{
		print(f3KAWho)
	}
	dat$plotMean( add=T )
	
	##
	#dat$plotMean(add=T)
}
dev.off()

#
png(sprintf('bioGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit = readRDS(fWho)
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat = readRDS(dWho)
	#
	par(mar=c(1,1,1,0))
	dat$plotQuan( function(B){B}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,dat$B0) )
	fit$plotQuan( function(B){B}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
}
dev.off()

#
png(sprintf('numGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit = readRDS(fWho)
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat = readRDS(dWho)
	#
	par(mar=c(1,1,1,0))
	dat$plotQuan( function(N){N}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,dat$N0) )
	fit$plotQuan( function(N){N}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
}
dev.off()

#
png(sprintf('avgGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit = readRDS(fWho)
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat = readRDS(dWho)
	#
	par(mar=c(1,1,1,0))
	dat$plotQuan( function(B,N){B/N}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]) ) #, ylim=c(0,dat$N0) )
	fit$plotQuan( function(B,N){B/N}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
}
dev.off()

#
png(sprintf('recGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit = readRDS(fWho)
	#
	f = function(t, Y, p=NULL){ der(t, Y, fit$lalpha, fit$lbeta, fit$gamma, fit$aS, fit$a0, fit$WW, fit$kappa, fit$catch, fit$B0) }
	fOut = dede(c(fit$N0, fit$B0), 1:TT, f, parms=NULL, method="lsode")#"radau") # 
	colnames(fOut) = c('time', 'N', 'B', 'Rec. Biomass', 'Net Growth')
	#
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat = readRDS(dWho)
	#
	d = function(t, Y, p=NULL){ der(t, Y, dat$lalpha, dat$lbeta, dat$gamma, dat$aS, dat$a0, dat$WW, dat$kappa, dat$catch, dat$B0) }
	dOut = dede(c(dat$N0, dat$B0), 1:TT, d, parms=NULL, method="lsode")#"radau") # 
	colnames(fOut) = c('time', 'N', 'B', 'Rec. Biomass', 'Net Growth')
	#
	par(mar=c(1,1,1,0))
	plot(dOut[,1], dOut[,4], 'l', lwd=3, ylim=c(0, max(dOut[,4], fOut[,4], na.rm=T)))
	lines(fOut[,1], fOut[,4], lwd=3, col='red')
	#dat$plotQuan( function(B,N){B/N}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]) ) #, ylim=c(0,dat$N0) )
	#fit$plotQuan( function(B,N){B/N}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
}
dev.off()

#
png(sprintf('groGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit = readRDS(fWho)
	#
	f = function(t, Y, p=NULL){ der(t, Y, fit$lalpha, fit$lbeta, fit$gamma, fit$aS, fit$a0, fit$WW, fit$kappa, fit$catch, fit$B0) }
	fOut = dede(c(fit$N0, fit$B0), 1:TT, f, parms=NULL, method="lsode")#"radau") # 
	colnames(fOut) = c('time', 'N', 'B', 'Rec. Biomass', 'Net Growth')
	#
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat = readRDS(dWho)
	#
	d = function(t, Y, p=NULL){ der(t, Y, dat$lalpha, dat$lbeta, dat$gamma, dat$aS, dat$a0, dat$WW, dat$kappa, dat$catch, dat$B0) }
	dOut = dede(c(dat$N0, dat$B0), 1:TT, d, parms=NULL, method="lsode")#"radau") # 
	colnames(fOut) = c('time', 'N', 'B', 'Rec. Biomass', 'Net Growth')
	#
	par(mar=c(1,1,1,0))
	plot(dOut[,1], dOut[,5], 'l', lwd=3, ylim=c(0, max(dOut[,4], fOut[,4], na.rm=T)))
	lines(fOut[,1], fOut[,5], lwd=3, col='red')
	#dat$plotQuan( function(B,N){B/N}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]) ) #, ylim=c(0,dat$N0) )
	#fit$plotQuan( function(B,N){B/N}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
}
dev.off()

#
png(sprintf('srrGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	#
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit  = readRDS(fWho)
	#
	wwFit = vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0) #fit$WW*(1-exp(-fit$kappa*fit$a0))
	FFit = FMsy(fit$M, fit$kappa, wwFit, fit$WW, exp(fit$lalpha), exp(fit$lbeta), fit$gamma)
	#
	fFit = function(x){wwFit*SRR(x, fit)-fit$M*x - fit$kappa*x + fit$kappa*fit$WW*NCar(x, 0, fit$M, fit$kappa, wwFit, fit$WW, exp(fit$lalpha), exp(fit$lbeta), fit$gamma) }
	mFit = stats::optimize(fFit, c(0, fit$B0), maximum=T)
	
	#
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat  = readRDS(dWho)
	#
	wwDat = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0) #dat$WW*(1-exp(-dat$kappa*dat$a0))
	FDat = FMsy(dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
	#
	fDat = function(x){wwDat*SRR(x, dat)-dat$M*x - dat$kappa*x + dat$kappa*dat$WW*NCar(x, 0, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma) }
	mDat = stats::optimize(fDat, c(0, dat$B0), maximum=T) 
	
	#
	fitFU = BBar(FFit, fit$M, fit$kappa, wwFit, fit$WW, exp(fit$lalpha), exp(fit$lbeta), fit$gamma)*FFit/mFit$objective
	datFU = BBar(FDat, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)*FDat/mDat$objective
	
	#
	par(mar=c(1,1,1,0))
	curve(fitFU*fDat(x), 0, max(dat$B0, fit$B0), n=10000, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), lwd=3, ylim=c(0, max(mFit$objective*fitFU, mDat$objective*datFU)))
	curve(fitFU*fFit(x), 0, fit$B0, col='red', add=T, lwd=3, n=10000)
	#netSRD = function(x){SRR(x, dat)-dat$M*x}
	#curve(netSRD(x), 0, dat$B0, lty=2, add=T)
	#netSRF = function(x){SRR(x, fit)-fit$M*x}
	#curve(netSRF(x), 0, fit$B0, lty=2, col='red', add=T)
	#curve(x*FDat, 0, dat$B0, add=T, col='red')
	#abline( v=BBar(FDat, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma) )
	##abline( v=dat$B0 )
	#curve(x*FFit, 0, fit$B0, add=T, col='red')
	#abline( v=BBar(FFit, fit$M, fit$kappa, wwFit, fit$WW, exp(fit$lalpha), exp(fit$lbeta), fit$gamma) )
	##dat$plotQuan( function(B,N){B/N}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]) ) #, ylim=c(0,dat$N0) )
	##fit$plotQuan( function(B,N){B/N}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
	#
	##
	#print( BBar(FFit, fit$M, fit$kappa, wwFit, fit$WW, exp(fit$lalpha), exp(fit$lbeta), fit$gamma)*FFit )
	#print( mFit$objective )
	#print( BBar(FFit, fit$M, fit$kappa, wwFit, fit$WW, exp(fit$lalpha), exp(fit$lbeta), fit$gamma)*FFit/mFit$objective )
	#print('')
	#print( BBar(FDat, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)*FDat )
	#print( mDat$objective )
	#print( BBar(FDat, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)*FDat/mDat$objective )
	#print('')
	#print('')
}
dev.off()

#
png(sprintf('rGrid%s.png', mod), width=2000, height=2000)
layout(lay)
for(i in 1:nrow(l)){ #nrow(out$ll)){
	#
	if( !i%in%tops ){next}
	#
	fWho = sprintf('%sfit_%s', place, rownames(l)[i])
	fit  = readRDS(fWho)
	#
	wwFit = vbGrow(fit$aS, fit$kappa, fit$WW, fit$a0) #fit$WW*(1-exp(-fit$kappa*fit$a0))
	fFit = function(x){SRR(x, fit)} #- fitt$kappa*x + fit$kappa*fit$WW*fit$N0
	mFit = stats::optimize(fFit, c(0, fit$B0), maximum=T)
	#
	dWho = sprintf('%sdatGen_%s', place, rownames(l)[i])
	dat  = readRDS(dWho)
	#
	wwDat = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0) #dat$WW*(1-exp(-dat$kappa*dat$a0))
	fDat = function(x){SRR(x, dat)} #- dat$kappa*x + dat$kappa*dat$WW*dat$N0}
	mDat = stats::optimize(fDat, c(0, dat$B0), maximum=T) 
	#
	par(mar=c(1,1,1,0))
	curve(fDat(x), 0, max(dat$B0, fit$B0), n=10000, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), lwd=3, ylim=c(0, max(mFit$objective, mDat$objective)))
	#curve(fFit(x), 0, fit$B0, col='red', add=T, lwd=3, n=10000)
	#FF = FMsy(dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
	#curve(x*FF, 0, dat$B0, add=T, col='red')
	#BB = BBar(FF, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
	#abline( v=BB )
	#abline( v=dat$B0 )
	#dat$plotQuan( function(B,N){B/N}, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]) ) #, ylim=c(0,dat$N0) )
	#fit$plotQuan( function(B,N){B/N}, add=T, col='red') #, main=sprintf("%1.2f, %1.2f", l[i,1], l[i,2]), ylim=c(0,fit$B0) ) #out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
}
dev.off()


##
#out = lhsForPlot(xiLim, zetaLim, N, 0, save=F)
##
#ox = order(out$ll[,1])
#lay = matrix(0, N, N) 
#for(i in 1:length(ox)){ lay[i,which(ox==i)]=i }
#lay = lay[seq(N,1),]

##
#png(sprintf('numGridA%sK%sW%s.png', a0, kappa, WW), width=2000, height=2000)
#layout(lay)
#for(i in 1:nrow(out$ll)){
#	par(mar=c(1,1,1,0))
#	out$mod[[i]]$plotQuan( function(N){N}, main=sprintf("%1.2f, %1.2f", out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$N0) )
#}
#dev.off()
#
##
#png(sprintf('avgGridA%sK%sW%s.png', a0, kappa, WW), width=2000, height=2000)
#layout(lay)
#for(i in 1:nrow(out$ll)){
#	par(mar=c(1,1,1,0))
#	out$mod[[i]]$plotQuan( function(B,N){B/N}, main=sprintf("%1.2f, %1.2f", out$ll[i,1], out$ll[i,2]) ) #, ylim=c(0,1) )
#}
#dev.off()
#
##
#png(sprintf('srrGridA%sK%sW%s.png', a0, kappa, WW), width=2000, height=2000)
#layout(lay)
#for(i in 1:nrow(out$ll)){
#	par(mar=c(1,1,1,0))
#	f = function(x){SRR(x, out$mod[[i]])-out$mod[[i]]$M*x}
#	curve(f(x), 0, out$mod[[i]]$B0, main=sprintf("%1.2f, %1.2f", out$ll[i,1], out$ll[i,2]) ) #, ylim=c(0,1) )
#}
#dev.off()





#dev.new()
#plot(out$ll)






	##record values
        #xiList[bin] = xi
        #zetaList[bin] = zs

##
#thresh = 0.02
#n = 28 #about 3 flushes all of the thialacia ranks
#n = 56

##
#if( F ){
##Only run this to start a design
#p = "./modsDDExpT45N150A1K10/" #'test'
#if(dir.exists(p)){ unlink(p, recursive=TRUE) }
#dir.create(p)
#
 #(xiLim, zetaLim, 
#}else{
#
##
#inPlace = "./modsDDExpT45N150A15/" 		#"./modsSchnuteHHardFlatT30N150WWideN84/"#"./modsSchnuteExpT30L3N150Wide/" #"./modsSchnuteHHardExpT45N150M0.1Wide/" #"./modsSchnuteHHardFlatT30
#outPlace = sprintf("./modsDDExpT45N150A15N%d/", n) 	#sprintf('./modsSchnuteHHardFlatT30N150WWideN%s/', 4*n) #sprintf('./modsSchnuteHHardFlatT30N150WWideAdapt%s/', thresh)
#
##
##MAIN
##
#
##
#datFiles = sprintf("%s%s", inPlace, list.files(path=inPlace, pattern=glob2rx("datGen*.rda")))
#fitFiles = sprintf("%s%s", inPlace, list.files(path=inPlace, pattern=glob2rx("fit*.rda")))
#
##
#oldMaxBatch = 0
##safely create place
#if( dir.exists(outPlace) ){
#        #
#        isGo = readline(sprintf("Overwrite %s Design?\nCtrl-C to Escape, Enter to Overwrite, or A to Append: ", outPlace))
#        #Append
#        if( toupper(isGo)=="A" ){
#                #Append
#                dir.create(outPlace, showWarnings=FALSE)
#                #update datFiles to come from the current modsPlace
#                datFiles = sprintf("%s%s", outPlace, list.files(path=outPlace, pattern=glob2rx("datGen*.rda")))
#                fitFiles = sprintf("%s%s", inPlace, list.files(path=inPlace, pattern=glob2rx("fit*.rda")))
#                oldMaxBatch = max( read.csv(sprintf("%s/appends.csv", outPlace))$batch )
#        #Overwrite      
#        }else{
#                #remove and recopy inPlace over to outPlace
#                unlink(outPlace, recursive=TRUE)
#                system(sprintf("cp -r %s %s", inPlace, outPlace))
#                #write.table(t(c('batch', 'modName')), sprintf("%s/appends.csv", outPlace), sep=',', row.names=F, col.names=F, append=F)
#                #dir.create(outPlace, showWarnings=FALSE)
#                #modsSchnuteHHardFlatT30N150WWide/ modsSchnuteHHardFlatT30N150WWideN28/")
#        }
##Make Fresh
#} else{
#        # 
#        system(sprintf("cp -r %s %s", inPlace, outPlace))
#        #write.table(t(c('batch', 'modName')), sprintf("%s/appends.csv", outPlace), sep=',', row.names=F, col.names=F, append=F)
#}
##dir.create(outPlace, showWarnings=FALSE) }
#
##
#xiList = c()
#zetaList = c()
#for(fitF in fitFiles){
#        #
#        fit = readRDS(fitF)
#        #check rsCOV
#        if(length(fit$rsCov)==0){ next }
#        #
#        xiList = c(xiList, fit$xi)
#        zetaList = c(zetaList, fit$zeta)
#}
#l = cbind(xiList, zetaList)
##vchull = l[chull(l),]
#
##
##png("adaptDesignSquare.png")
#if( !is.null(l) ){ plot(l[,1], l[,2]) }
#
##
#ii = 1
#ms = c()
##
#for(i in 1:n){
#	#
#	mOut = maximinAddHullScale(1, xlim, ylim, XStart=l)
#        m = mOut$X
#	m01 = mOut$Z
#        r01 = mOut$md
#        #
#        X = rbind(l, m)
#        ms = c(ms, r01)
#        #
#        print(sprintf("(%3s, %3s): %s", m01[1], m01[2], r01))
#        #       
#        lx = move01tolim(m01[1]-r01, xlim)
#        ux = move01tolim(m01[1]+r01, xlim)
#        ly = move01tolim(m01[2]-r01, ylim)
#        uy = move01tolim(m01[2]+r01, ylim)	
#	#
#        xB = c(max(lx, xlim[1]), min(ux, xlim[2]))
#        yB = c(max(ly, ylim[1]), min(uy, ylim[2]))
#        #
#        ll = lhsMake(xB, yB, 5, ii+oldMaxBatch, save=outPlace)
#        l = rbind(l, ll)
#
#        #
#        text(m, col=i, label=sprintf("%s",i))
#        points(ll, pch=19, col=i) #cols[i])
#        polygon(l[chull(l),], border=i)
#        segments(xB[1], yB[1], xB[2], yB[1], col=i, lty=3)
#        segments(xB[1], yB[2], xB[2], yB[2], col=i, lty=3)
#        segments(xB[1], yB[1], xB[1], yB[2], col=i, lty=3)
#        segments(xB[2], yB[1], xB[2], yB[2], col=i, lty=3)
#
#        #
#	#print(nrow(ll))
#        if( nrow(ll)>0 ){ ii=ii+1 }
#}
#
#}









##NOTE: Only Temporary for testing
#xiList = c()
#zetaList = c()
#for(fitF in datFiles){
#        #
#        fit = readRDS(fitF)
#        #check rsCOV
#        #if(length(fit$rsCov)==0){ next }
#        #
#        xiList = c(xiList, fit$xi)
#        zetaList = c(zetaList, fit$zeta)
#}
#l = cbind(xiList, zetaList)





##
#xiLim = c(0.5, 3.75)
#zetaLim = c(0.15, 0.7)



##
#Fmsy  = FMsy(M, kappa, ww, WW, alpha, beta, gamma)
#BStar = BBar(Fmsy, M, kappa, ww, WW, alpha, beta, gamma)
#BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)



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



