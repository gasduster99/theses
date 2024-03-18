rm(list=ls())

#
library(plgp)
library(crch)
library(gMOIP)
library(Ryacas)
library(pracma)
library(stringr)
library(latex2exp)
library(rootSolve)
library(RColorBrewer)

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
getZetaBH = function(x, M, k, W, a0){
        #
        w = W*(1-exp(-k*a0))
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

#
#MAIN
#

#
a0 = 0.1 	#15     #15  #7.5 #15  #1
M  = 0.2
kappa = 15 	#0.1 #0.2 #0.2 #10
WW = 1
ww = WW*(1-exp(-kappa*a0))

#
B0 = 10000

#just initiating to give alpha, beta, gamma some reasonable values
xi = 1
zeta = 0.4
#
inv = invert(zeta, xi, B0, M, kappa, ww, WW)
alpha = inv$alpha
beta  = inv$beta
gamma = inv$gamma

#
xis = seq(0.05, 4, 0.1)
#
aI = 80 #15#7 #15   #10 #
kI = 80 #15#15   #10 #
#mI = 150#150  #150#
#wI = 15
as = logseq(0.01, 10, n=aI) #seq(0.01, 15, length.out=aI) #0.5)
ks = logseq(0.1, 2, n=kI) #seq(0.01, 15, length.out=kI) #0.5)
#as = seq(0.01, 10, len=aI)
#ks = seq(0.01, 10, len=kI)
#Ms = logseq(0.01, 1, n=mI)  #seq(0.05, 1, length.out=mI)  #0.1)
#Ws = logseq(0.01, 15, n=wI)
#
cols = brewer.pal(5, 'Set1')
kCol = c(255, 0, 0) #col2rgb(cols[1])
aCol = c(0, 0, 255) #c(0, 255, 0) #col2rgb(cols[2])
wwCol = colorRamp(c(rgb(aCol[1]/255, aCol[2]/255, aCol[3]/255), rgb(kCol[1]/255, kCol[2]/255, kCol[3]/255)))
#mCol = c(0, 0, 255) #col2rgb(cols[3])
#wCol = mCol
#
#xiM = matrix(rep(xis, length(xis)), length(xis), length(as), byrow=T)
#bhL = c() #matrix(NA, nrow=length(as), ncol=length(xis))
#plot(xis, getZetaBH(xis, mean(Ms), mean(ks), WW, mean(as)), ylim=c(0.2, 0.7), type='l')
#png("rpSpace.png", width=700, height=700)
pdf("rpSpaceww.pdf")
plot(xis, 1/(xis+2), col='white', type='l', lwd=2, ylim=c(0.1, 0.6), 
	xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'),  
	main="Space of BH Reference Points"
)
#
for(ai in 1:aI){ #sample(1:aI)){	
	for(ki in 1:kI){ #sample(1:kI)){
		#for(mi in sample(1:mI)){
		#for(wi in sample(1:wI))
			#
			a = as[ai]
			k = ks[ki]
			m = 0.2 #Ms[mi]
			ww = WW*(1-exp(-k*a))
			#m = 0.2
			#w = Ws[wi]
			#
			bhLine <- tryCatch({
				bhLine = getZetaBH(xis, m, k, WW, a)
			}, error=function(err){
				bhLine = rep(NA, length(xis))
			})
			#bhL = cbind(bhL, bhLine)
			mixCol = rgb(wwCol(ww)[1]/255, wwCol(ww)[2]/255, wwCol(ww)[3]/255, alpha=0.1)
			#mixCol = rgb(
			#	(ki/kI)*kCol[1], 
			#	0,
			#	(ai/aI)*aCol[3],
			#	#0, 
			#	#(mi/mI)*mCol[3],
			#	#(wi/wI)*wCol[3],
			#	alpha=3*255/(aI+kI), #+mI), 
			#	maxColorValue=255
			#)
			##print(mixCol)
			lines(xis, bhLine, col=mixCol, lwd=5)
			#print(c(a,k))
		#}
	}
}
#
lines(xis, 1/(xis+2), lwd=2)
legend("topright", 
	legend = c(
		TeX("$w(a_s)$ High"), TeX("$w(a_s)$  Low")
		#, TeX("$M$ High")
	), 
	col = c(
		rgb(255,0,0,alpha=255, maxColorValue=255), rgb(0,0,255,alpha=255, maxColorValue=255) 
		#,rgb(0,0,255,alpha=100, maxColorValue=255)
	), 
	lwd = 5
)
dev.off()



##
#matlines(t(xiM), t(bhL))
#((ai/aI)*aCol[1] + (ki/kI)*kCol[1] + (mi/mI)*mCol[1])/3, 
				#((ai/aI)*aCol[2] + (ki/kI)*kCol[2] + (mi/mI)*mCol[2])/3, 
				#((ai/aI)*aCol[3] + (ki/kI)*kCol[3] + (mi/mI)*mCol[3])/3, 









#OLD

##DIFFERENTIAL EQUATION 
#
##
#der = function(t, Y, lalpha, lbeta, gamma, a0, WW, kappa, catch, B0){
#        #linearly interpolate catches
#        ft = floor(t)
#        q  = (t-ft)
#        Cl = catch[ft]
#        Cu = catch[ft+1]
#        C = q*Cu + (1-q)*Cl
#        if(q==0){ C=Cl }
#        #
#        N = Y[1]
#        B = Y[2]
#        #
#        if( (t-a0)<1){
#                Blag = B0
#        }else{
#                Blag = lagvalue(t-a0)[2]
#        }
#        #
#        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
#        alpha = exp(lalpha)
#        beta = exp(lbeta)
#        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
#        #
#	ww = WW*(1-exp(-kappa*a0))
#	FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
#        #
#        out = c(N=NA, B=NA)
#        out[1] = R - (M+FF)*N
#        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B 
#        #
#        return( list(out) )
#}
#
##
#SRR = function(B, mod){
#	#
#	alpha = exp(mod$lalpha)
#        beta = exp(mod$lbeta)
#	gamma = mod$gamma
#	#
#        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
#	#
#	return(R)
#}
#SRR = Vectorize(SRR, 'B')
#
##DESIGN STUFF
#
##
#lhsForPlot = function(xiLim, zetaLim, n, batch, save=F){
#        #somehow define place when save!=F
#        #deal with append behaviour	
#	
#	#LHS boundaries
#	xiMin = xiLim[1]        #0.25 #3.25 #0.5
#	xiMax = xiLim[2]        #3.75 #3.5 
#	zetaMin = zetaLim[1]    #0.15 #0.32 #0.15 
#	zetaMax = zetaLim[2]    #0.7  #0.43 #0.7 
#	#zeta range should not be too small
#	
#	#xi and zeta Bins (bin defined by left edge right egde not used) [n+1 used to give n left 
#	xiL = seq(xiMin, xiMax, length.out=n+1)
#	zetaL = seq(zetaMin, zetaMax, length.out=n+1)
#	minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
#	binTrk = ceiling(abs(log10(minDiff)))
#	#
#	mod = vector(mode='list', length=n)
#	Llist = rep(T, n) #a list to track which zeta bin still needs a sample
#	zetaList = rep(NA, n) #a list of the sampled zeta value in each bin
#	xiList = rep(NA, n)   #a list of the sampled xi value in each bin	
#
#	#details of how to sample the curvy bit for comuting sd
#	scope = 2 #how far to sample the inverse function for curvature
#	scopeN = 30 #how many samples to take of the inverse function curvature
#	scopeRun = seq(-scope, scope, length.out=scopeN)
#	#start at fudge (hessian sd) slowing increase fudge by fudgeFactor so that we eventually s
#	fudge = 1
#	giveUp = 100 		#100 #10
#	fudgeFactor = 0.1 	#0.1 #0.01
#	#start df at dfStart with thin tails (to explore curvature), decrement by dfFactor to tran
#	dfStart = 100
#	df = dfStart
#	dfFactor = 1 #dfFactor, dfStart, and fudgeFactor are balanced to first decrease df, and th
#	
#	#
#	xii = 1
#	xiShuffle = sample(xiL[-(n+1)])
#	#
#	while( any(Llist) ){
#		#for each xi column (n of them)
#		xiLWho = which(xiL==xiShuffle[xii])
#		xi = runif(1, xiShuffle[xii], xiL[which(xiL==xiShuffle[xii])+1])
#		ff = xi*M
#		
#		#should return NA if something fails so that minZ works as planned
#		z = function(g){ 
#			a = getAlphaFmsy(ff, M, kappa, ww, WW, beta, g)
#			getZeta(ff, M, kappa, ww, WW, a, beta, g) 
#		}
#		z = Vectorize(z, "g")
#		minZ = suppressWarnings(stats::optimize(z, c(-100, 100))$minimum)
#		f = function(x){grad(z, x)}#pdf analogy    
#		fv = Vectorize(f, "x")
#		peakZ = suppressWarnings(stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))$minimum)
#		varZ = tryCatch({ 1/optimHess(peakZ, function(x){-log(fv(x))})
#		}, error=function(err){ varZ=-1 })
#		if(varZ<0){ 
#			#
#			fudge = fudge+fudgeFactor #increase fudge factor to make progress toward giving up
#                        if( fudge>giveUp ){ break }
#			#
#			next 
#		}
#	        #
#		gs = rtt(1, peakZ, sqrt(varZ)*fudge, df, left=minZ)
#	
#		#propose sample of zeta
#		zs = z(gs)
#		#reject: catch the case where zs is propsed outside of z range
#		if(zs<zetaMin | zs>zetaMax){
#			# 
#			fudge = fudge+fudgeFactor #increase fudge factor to make progress toward giving up
#                        if( fudge>giveUp ){ break }
#			#
#			next 
#		}#NOTE: if slow, I could probably use the invert function to include bounds in truncation above	
#		
#		#
#		as = getAlphaFmsy(ff, M, kappa, ww, WW, beta, gs)
#		bs = getBeta(B0, M, kappa, ww, WW, as, gs)
#	
#		#
#		fStar = FMsy(M, kappa, ww, WW, as, bs, gs)
#		#if( length(fStar)==0 ){
#		#	#print('Skip Fmsy')
#		#	next
#		#}
#		
#		#find left edge of bin sampled
#		bin = max(which(zetaL<=zs))
#		#if sampled bin not yet sampled, then save
#		if( Llist[bin] ){
#		        #record values
#		        xiList[bin] = xi
#		        zetaList[bin] = zs
#			
#			#
#			dat = ddModel$new( derivs=der,
#				N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){	#virgin numbers
#				        #
#				        alpha = exp(lalpha)
#				        beta  = exp(lbeta)
#					ww = WW*(1-exp(-kappa*a0))
#				        #
#					BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma) 
#				        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
#					
#				},
#				B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){	#virgin biomass
#				        #
#				        alpha = exp(lalpha)
#				        beta = exp(lbeta)
#					ww = WW*(1-exp(-kappa*a0))
#				        #
#					BBar(0, M, kappa, ww, WW, alpha, beta, gamma) 	
#				},
#				time=1:TT, catch=FtFmsy, a0=a0, M=M, WW=WW, kappa=kappa,#constants
#				lalpha=log(as), lbeta=log(bs), gamma=gs, 		#parameters
#				lq=log(0.00049), lsdo=log(0.01160256),	        	#nuisance parameters
#				xi=xi, zeta=zs
#			)
#			dat$iterate(odeMethod)
#			mod[[bin]] = dat
#			#dat$plotQuan(main=sprintf("%s, %s", xi, zs))
#			#writeLines(sprintf("\n\n%s, %s", xi, zs))
#			#dat$printSelf()
#			
#			#       
#                        if( save!=F ){
#                                if( save==T ){
#                                        #create place
#                                        place = date()
#					save = place
#					dir.create(place, showWarnings=FALSE)
#                                }else{ place=save }
#                                #NOTE: catelogue by bin lables (left edges of xiL, zetaL)s  xiShuffle[xii]
#                                #datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiL[bin], binTrk), round(zetaL[bin], binTrk)))
#                                datName = sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiShuffle[xii], binTrk), round(zetaL[bin], binTrk))
#                                dat$save(datName)
#				#
#				apPlace = sprintf("%s/appends.csv", place)
#				if( !file.exists(apPlace) ){ 
#					write.table(t(c('batch', 'modName')), apPlace, sep=',', row.names=F, col.names=F, append=F)
#				}
#                                write.table(t(c(batch, datName)), apPlace, sep=',', row.names=F, col.names=F, append=T)
#                        }
#			
#			#knobs for making the design
#			Llist[bin] = F  #turn off bin
#			xii = xii+1     #increment xi counter
#			df = 100        #reset df to thin tails
#			fudge = 1       #reset fudge to thin tails
#		}else{ # catch case where hyperbola gap requires heavy tails
#			df = max(1, df-1)  #decrement df to broadend search to get the center bald spot
#			fudge = fudge+fudgeFactor #increase fudge factor to get a wider sd
#			if( fudge>giveUp ){ break } #give up (it mostly just misses one or two bins so its usually fine)
#			#writeLines(sprintf('\t%s', fudge))
#		}
#	}
#	#
#	out = list(
#		ll = cbind(xiList[!is.na(xiList)], zetaList[!is.na(zetaList)]),
#		mod = mod[!is.na(xiList)]
#	)
#	#
#	return( out )
#}
#
##
#move01tolim = function(coord01, lim){
#        coord01*(lim[2]-lim[1])+lim[1]
#}
#
##
#movelimto01 = function(coordlim, lim){
#        (coordlim-lim[1])/(lim[2]-lim[1])
#}
#
##
#maximinAddHullScale = function(n, xlim, ylim, TT=1e5, XStart=NULL) {
#        #n : number of points to draw
#        #m : dimension of space to draw from
#        #
#        #value :  only returns new additions
#
#        #
#        ZStart = cbind( movelimto01(XStart[,1], xlim), movelimto01(XStart[,2], ylim) )
#        vchull = ZStart[chull(ZStart),]
#        #plot(ZStart)
#        #polygon(vchull)
#
#        #
#        Z = matrix(NA, nrow=n, ncol=2)
#        Z[,1] = rep(vchull[1,1], n) ## start everything in a place that is guaranteed to be inside the hull
#        Z[,2] = rep(vchull[1,2], n)
#        #
#        d = distance(Z)
#        d = d[upper.tri(d)]
#        #
#        md = min(d)
#        #handle initiated design
#        if(!is.null(ZStart)){
#                md2 = min(distance(Z, ZStart))
#                if(md2 < md){ md=md2 }
#        }
#        #
#        propZ = cbind(runif(TT), runif(TT))
#        isIn = inHull(propZ, vchull)>=0
#        inZ = propZ[isIn,]
#        #
#        tList = 1:sum(isIn)
#        for(t in tList){
#                #propse a row
#                row = sample(1:n, 1)
#                zold = Z[row,]
#                #assume we accept
#                Z[row,] = inZ[t,] ## random new row 
#                d = distance(Z) #new distance
#                d = d[upper.tri(d)]
#                mdprime = min(d)
#                #handle initial design (make min dist smaller if neccesary with DStart included)
#                if(!is.null(ZStart)){
#                        mdprime2 = min(distance(Z, ZStart))
#                        if(mdprime2 < mdprime) mdprime=mdprime2
#                }
#                #
#                if( mdprime > md ){ md=mdprime ## accept 
#                } else{ Z[row,]=zold } ## reject
#
#        }
#        #points(Z, pch=19)
#        #addCircle(Z[,1], Z[,2], sqrt(md))
#
#        #
#        X = cbind( move01tolim(Z[,1], xlim), move01tolim(Z[,2], ylim) )
#        out = list(X=X, Z=Z, md=sqrt(md))
#
#        #
#        return( out )
#}
#
##
#addCircle = function(centerx, centery, radius, length=200){
#        # prepare "circle data"
#        theta = seq(0, 2 * pi, length=length) # angles for drawing points around the circle
#
#        # draw the circle
#        lines(x = radius * cos(theta) + centerx, y = radius * sin(theta) + centery)
#}
#
#
##
##HEAD
##
#
##DD MODEL STUFF
#
##
#odeMethod = "lsode"
#
##
#a0 = 0.1 	#15     #15  #7.5 #15  #1
#M  = 0.2
#kappa = 15 	#0.1 #0.2 #0.2 #10
#WW = 1
#ww = WW*(1-exp(-kappa*a0))
#
##
#B0 = 10000
#
##just initiating to give alpha, beta, gamma some reasonable values
#xi = 1
#zeta = 0.4
##
#inv = invert(zeta, xi, B0, M, kappa, ww, WW)
#alpha = inv$alpha
#beta  = inv$beta
#gamma = inv$gamma
#
##
#TT = 45
#FtFmsy = rep(1, TT) #make faux catch
#
##DESIGN STUFF
#
##
#thresh = 0.02
#n = 28 #about 3 flushes all of the thialacia ranks
##n = 56
#
##
#xlim = c(0.25, 3.75)
#ylim = c(0.15, 0.7) #0.6) #
#
##NOTE:high zeta and high xi lets the initial drop get below Bmsy
#N = 20
##set.seed(111)
##
#out = lhsForPlot(xlim, ylim, N, 0, save=F)
#N = nrow(out$ll)
##
#ox = order(out$ll[,1])
#lay = matrix(0, N, N) 
#for(i in 1:length(ox)){ lay[i,which(ox==i)]=i }
#lay = lay[seq(N,1),]
#
##
#png(sprintf('bioGridA%sK%sW%s.png', a0, kappa, WW), width=2000, height=2000)
#layout(lay)
#for(i in 1:nrow(out$ll)){
#	par(mar=c(1,1,1,0))
#	out$mod[[i]]$plotQuan( function(B){B}, main=sprintf("%1.2f, %1.2f", out$ll[i,1], out$ll[i,2]), ylim=c(0,out$mod[[i]]$B0) )
#}
#dev.off()
#
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



