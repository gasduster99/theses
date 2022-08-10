rm(list=ls())

#
library(GA)
library(crch)
library(plgp)
library(gMOIP)
library(deldir)
library(pracma)
library(numDeriv)
library(rootSolve)
library(RColorBrewer)
#
source('prodClass0.1.1.r')

#make convex hull
#https://www.rdocumentation.org/packages/grDevices/versions/3.6.2/topics/chull
#maximin accept reject scheme, rejecting outisdie of the hull
#https://search.r-project.org/CRAN/refmans/gMOIP/html/inHull.html
#https://askubuntu.com/questions/913493/gsl-config-not-found

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
lhsMake = function(xiLim, zetaLim, n, batch, save=F){
	#somehow define place when save!=F
	#deal with append behaviour
	
	#
	TT=30
	FtFmsy = rep(1, TT) #make faux catch
	
	#
	B0 = 10000
	M  = 0.2
	odeMethod = "lsode"
	
	#number of samples in design
	#n = 150 #500 maxed out method
	
	#LHS boundaries
	xiMin = xiLim[1]	#0.25 #3.25 #0.5
	xiMax = xiLim[2]	#3.75 #3.5 
	zetaMin = zetaLim[1]	#0.15 #0.32 #0.15 
	zetaMax = zetaLim[2]	#0.7  #0.43 #0.7 
	#zeta range should not be too small
	
	#xi and zeta Bins (bin defined by left edge right egde not used) [n+1 used to give n left edges]
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
	#start at fudge (hessian sd) slowing increase fudge by fudgeFactor so that we eventually search flat region, but still allow for repetition at calculated sd 
	fudge = 1
	giveUp = 10 #100
	fudgeFactor = 0.1 #0.01
	#start df at dfStart with thin tails (to explore curvature), decrement by dfFactor to transition to explore flat regiem 
	dfStart = 100
	df = dfStart
	dfFactor = 1 #dfFactor, dfStart, and fudgeFactor are balanced to first decrease df, and then crank on sd
	
	#
	xii = 1
	xiShuffle = sample(xiL[-(n+1)])
	#
	while( any(Llist) ){
		#for each xi column (n of them)
		xi = runif(1, xiShuffle[xii], xiL[which(xiL==xiShuffle[xii])+1])
		ff = xi*M
		
		#
		minZ = suppressWarnings(stats::optimize(function(x){z(x, ff, M)}, c(-100, 100))$minimum)
		f = function(x){grad(z, x, ff=ff, M=M)}#pdf analogy    
		fv = Vectorize(f, "x")
		peakZ = stats::optimize(function(x){-log(fv(x))}, c(minZ, 100))$minimum
		varZ = 1/optimHess(peakZ, function(x){-log(fv(x))})
		if(varZ<0){ next }
		#
		gs = rtt(1, peakZ, sqrt(varZ)*fudge, df, left=minZ)
		
		#propose sample of zeta
		zs = z(gs, ff, M)
		#reject: catch the case where zs is propsed outside of z range
		if(zs<zetaMin | zs>zetaMax){ next } 
		
		#NOTE: compute auxilliary quantities
		#NOTE: make sure that reverse a,b,g->xi,zeta and xi,zeta->a,b,g calcs make sense; else reject	
		as = a(gs, ff, M)
		bs = b(as, gs, M, B0)
		
		#
		fStar = FMsy(as, gs, M)
		if( length(fStar)==0 ){ 
			#print('Skip Fmsy')
			next
		}
			
		#find left edge of bin sampled
		bin = max(which(zetaL<=zs))
		#if sampled bin not yet sampled, then save
		if( Llist[bin] ){
			#record values
			xiList[bin] = xi 
			zetaList[bin] = zs
			#
			#NOTE: generate a faux data object and catelogue design elements in a model object in the design directory
			#make Schunte production model
	               	datGen = prodModel$new(
	               	        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
	               	        time=1:30, catch=FtFmsy, M=M,                   #constants
	               	        alpha=as, beta=bs, gamma=gs,                   #parameters
	               	        lalpha=log(as), lbeta=log(bs),                 #reparameterize
	               	        lq=log(0.00049), lsdo=log(0.01160256), #log(0.1160256) #nuisance paramet
	               	        xi=xi, zeta=zs                                #other incidentals to carr
	               	)
	               	datGen$iterate(odeMethod)
			#datGen$plotMean()
			
			#	
			if( save!=F ){
				if( save==T ){ 
					#create place
					place=date(); dir.create(place, showWarnings=FALSE)
				}else{ place=save }
				#NOTE: catelogue by bin lables (left edges of xiL, zetaL)s  xiShuffle[xii]
				#datGen$save(sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiL[bin], binTrk), round(zetaL[bin], binTrk)))
				datName = sprintf('%s/datGen_xi%s_zeta%s.rda', place, round(xiShuffle[xii], binTrk), round(zetaL[bin], binTrk))
				datGen$save(datName)	
				write.table(t(c(batch, datName)), sprintf("%s/appends.csv", outPlace), sep=',', row.names=F, col.names=F, append=T)
			}

			#knobs for making the design
			Llist[bin] = F 	#turn off bin
			xii = xii+1 	#increment xi counter
			df = 100    	#reset df to thin tails
			fudge = 1   	#reset fudge to thin tails
		}else{ # catch case where hyperbola gap requires heavy tails
			df = max(1, df-1)  #decrement df to broadend search to get the center bald spot
			fudge = fudge+fudgeFactor #increase fudge factor to get a wider sd
			if( fudge>giveUp ){ break } #give up (it mostly just misses one of two bins so its usually fine)
		}	
	}
	#
	#return( list(xiList=xiList, zetaList=zetaList) )
	return( cbind(xiList[!is.na(xiList)], zetaList[!is.na(zetaList)]) )
}

#
maximinAdd = function(n, xlim, ylim, TT=5e5, XStart=NULL) {
        #n : number of points to draw
        #m : dimension of space to draw from
        #
        #value :  only retunrs new additions
        
	#
        X = matrix(NA, nrow=n, ncol=2)
	X[,1] = runif(n, xlim[1], xlim[2])   ## random new row
	X[,2] = runif(n, ylim[1], ylim[2])
	d = distance(X)
        d = d[upper.tri(d)]
        #
        md = min(d)
        #handle initiated design
        if(!is.null(XStart)){
                md2 = min(distance(X, XStart))
                if(md2 < md){ md=md2 }
        }
        for(t in 1:TT){
                #propse a row
                row = sample(1:n, 1)
                xold = X[row,]
                #assume we accept
                #X[row,] = runif(m)
		X[row,1] = runif(1, xlim[1], xlim[2])   ## random new row 
                X[row,2] = runif(1, ylim[1], ylim[2])
                d = distance(X) #new distance
                d = d[upper.tri(d)]
                mdprime = min(d)
                #handle initial design (make min dist smaller if neccesary with DStart included)
                if(!is.null(XStart)){
                        mdprime2 = min(distance(X, XStart))
                        if(mdprime2 < mdprime) mdprime=mdprime2
                }
                #
                if( mdprime > md ){ md=mdprime ## accept 
                } else{ X[row,]=xold } ## reject
        }
        #return only the augmenting points
        return(X)
}

#
maximinAddHull = function(n, vchull, TT=5e5, XStart=NULL) {
        #n : number of points to draw
        #m : dimension of space to draw from
        #
        #value :  only returns new additions
        
	#
	xlim = c(min(vchull[,1]), max(vchull[,1]))
	ylim = c(min(vchull[,2]), max(vchull[,2]))
	#
        X = matrix(NA, nrow=n, ncol=2)
	X[,1] = rep(mean(xlim), n) ## start everything in a place that is guaranteed to be inside the hull
	X[,2] = rep(mean(ylim), n) 
	d = distance(X)
        d = d[upper.tri(d)]
        #
        md = min(d)
        #handle initiated design
        if(!is.null(XStart)){
                md2 = min(distance(X, XStart))
                if(md2 < md){ md=md2 }
        }
	propX = cbind(runif(TT, xlim[1], xlim[2]), runif(TT, ylim[1], ylim[2]))
	isIn = inHull(propX, vchull)>=0
	inX = propX[isIn,]
	tList = 1:sum(isIn)
        for(t in tList){
                #propse a row
                row = sample(1:n, 1)
                xold = X[row,]
                #assume we accept
                X[row,] = inX[t,] ## random new row 
                d = distance(X) #new distance
                d = d[upper.tri(d)]
                mdprime = min(d)
                #handle initial design (make min dist smaller if neccesary with DStart included)
                if(!is.null(XStart)){
                        mdprime2 = min(distance(X, XStart))
                        if(mdprime2 < mdprime) mdprime=mdprime2
                } 
		#
		if( mdprime > md ){ md=mdprime ## accept 
                } else{ X[row,]=xold } ## reject
        	
	}
        #return only the augmenting points
        return(X)
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
#MAIN
#

#
thresh = 0.02 
r01 = 1
n = 28 #56 #about 3 flushes all of the thialacia ranks
#r = 1
#r1 = r
#t = 1
#
xlim = c(0.25, 3.75)
ylim = c(0.15, 0.7)

#
inPlace = "./modsSchnuteExpT30L3N150Wide/" #"./modsSchnuteHHardExpT45N150M0.1Wide/" #"./modsSchnuteHHardFlatT30N150WWideN56/" #'./modsSchnuteHHardFlatT30N150WWide/' #Extra/'
outPlace = sprintf('./modsSchnuteHHardExpT30L3N150WWideN%s/', n) #3*n) #sprintf('./modsSchnuteHHardFlatT30N150WWideAdapt%s/', thresh)
#
datFiles = sprintf("%s%s", inPlace, list.files(path=inPlace, pattern=glob2rx("datGen*.rda")))
fitFiles = sprintf("%s%s", inPlace, list.files(path=inPlace, pattern=glob2rx("fit*.rda")))

#start append file
#remove overwrite options for now
	#later read append file to roll back appends
#add in a line to add to the append file

#
oldMaxBatch = 0
#safely create place
if( dir.exists(outPlace) ){
        #
        isGo = readline(sprintf("Overwrite %s Design?\nCtrl-C to Escape, Enter to Overwrite, or A to Append: ", outPlace))
        #Append
        if( toupper(isGo)=="A" ){
                #Append
                dir.create(outPlace, showWarnings=FALSE)
		#update datFiles to come from the current modsPlace
		datFiles = sprintf("%s%s", outPlace, list.files(path=outPlace, pattern=glob2rx("datGen*.rda")))
		fitFiles = sprintf("%s%s", inPlace, list.files(path=inPlace, pattern=glob2rx("fit*.rda")))
		oldMaxBatch = max( read.csv(sprintf("%s/appends.csv", outPlace))$batch )
	#Overwrite  	
	}else{
                #remove and recopy inPlace over to outPlace
                unlink(outPlace, recursive=TRUE)
		system(sprintf("cp -r %s %s", inPlace, outPlace))
		write.table(t(c('batch', 'modName')), sprintf("%s/appends.csv", outPlace), sep=',', row.names=F, col.names=F, append=F)
                
		#dir.create(outPlace, showWarnings=FALSE)
		#modsSchnuteHHardFlatT30N150WWide/ modsSchnuteHHardFlatT30N150WWideN28/")
        }
#Make Fresh
} else{
	# 
	system(sprintf("cp -r %s %s", inPlace, outPlace)) 
	write.table(t(c('batch', 'modName')), sprintf("%s/appends.csv", outPlace), sep=',', row.names=F, col.names=F, append=F)
}
#dir.create(outPlace, showWarnings=FALSE) }

#
xiList = c()
zetaList = c()
for(fitF in fitFiles){
        #
        fit = readRDS(fitF)
        #check rsCOV
	if(length(fit$rsCov)==0){ next }
	#
        xiList = c(xiList, fit$xi)
        zetaList = c(zetaList, fit$zeta)
}
l = cbind(xiList, zetaList)
#vchull = l[chull(l),]

#
#png("adaptDesignSquare.png")
plot(l)

#
i = 1
ms = c()
#
for(i in 1:n){
#while(r01>thresh){
	#
	#m = maximinAddHull(1, vchull, XStart=l) #rbind(l, boxGrid))	
	mOut = maximinAddHullScale(1, xlim, ylim, XStart=l)
	m = mOut$X
	m01 = mOut$Z
	r01 = mOut$md
	#
	X = rbind(l, m)
	ms = c(ms, r01)
	#
	print(sprintf("(%3s, %3s): %s", m01[1], m01[2], r01))
	#	
	lx = move01tolim(m01[1]-r01, xlim)  
	ux = move01tolim(m01[1]+r01, xlim)  
	ly = move01tolim(m01[2]-r01, ylim)  
	uy = move01tolim(m01[2]+r01, ylim)  
	#
	xB = c(max(lx, xlim[1]), min(ux, xlim[2]))
	yB = c(max(ly, ylim[1]), min(uy, ylim[2]))
	#
	ll = lhsMake(xB, yB, 5, i+oldMaxBatch, save=outPlace)
	l = rbind(l, ll)
	
	#
	text(m, col=i, label=sprintf("%s",i))
	points(ll, pch=19, col=i) #cols[i])
	polygon(l[chull(l),], border=i)
	segments(xB[1], yB[1], xB[2], yB[1], col=i, lty=3)
	segments(xB[1], yB[2], xB[2], yB[2], col=i, lty=3)
	segments(xB[1], yB[1], xB[1], yB[2], col=i, lty=3)
	segments(xB[2], yB[1], xB[2], yB[2], col=i, lty=3)
	
	#
	i = i+1
	#t = r1-r
	#r1 = r
}
#dev.off()









#
#l = lhsMake(c(0.25, 3.75), c(0.15, 0.7), 20)
#l = matrix(unlist(l), ncol=2)


#xr = min(xDD)
#yr = min(yDD)
#l = lhsMake(c(0.25, 3.75), c(0.15, 0.7), 5)

##
#nn = 10
#boxGrid = cbind(seq(0.25, 3.75, length.out=nn*2), rep(0.15, nn*2))
#boxGrid = rbind(boxGrid, cbind(seq(0.25, 3.75, length.out=nn*2), rep(0.7, nn*2)) )
#boxGrid = rbind(boxGrid, cbind(rep(0.25, nn*0.5), seq(0.15, 0.7, length.out=nn*0.5)) )
#boxGrid = rbind(boxGrid, cbind(rep(3.75, nn*0.5), seq(0.15, 0.7, length.out=nn*0.5)) )
##
###
##tri = tri.mesh(l$x, l$z)
##
#nn = 10
#boxGrid = cbind(seq(0.25, 3.75, length.out=nn*2), rep(0.15, nn*2))
#boxGrid = rbind(boxGrid, cbind(seq(0.25, 3.75, length.out=nn*2), rep(0.7, nn*2)) )
#boxGrid = rbind(boxGrid, cbind(rep(0.25, nn*0.5), seq(0.15, 0.7, length.out=nn*0.5)) )
#boxGrid = rbind(boxGrid, cbind(rep(3.75, nn*0.5), seq(0.15, 0.7, length.out=nn*0.5)) )
###
##tess = deldir(c(l$x, boxGrid[,1]), c(l$z, boxGrid[,2]))
##tiles = tile.list(tess, minEdgeLength=min(distance(l$x, l$z)))















##
#maximin = function(n, xlim, ylim, TT=10^5, Xorig=NULL){ 
#	#
#	X <- matrix(runif(n*2), ncol=2) ## initial design 
#	X[,1] = X[,1]*(xlim[2]-xlim[1]) + xlim[1]
#	X[,2] = X[,2]*(ylim[2]-ylim[1]) + ylim[1]
#	d <- distance(X) 
#	d <- d[upper.tri(d)] 
#	md <- min(d) 
#	if( !is.null(Xorig) ){ 		## new code 
#		md2 <- min(distance(X, Xorig)) 
#		if(md2 < md) md <- md2
#	}
#	for( t in 1:TT ){ 
#		row <- sample(1:n, 1) 
#		xold <- X[row,] 	## random row selection 
#		X[row,1] = runif(1, xlim[1], xlim[2]) 	## random new row 
#		X[row,2] = runif(1, ylim[1], ylim[2])
#		d <- distance(X) 
#		d <- d[upper.tri(d)] 
#		mdprime <- min(d) 
#		if(!is.null(Xorig)){ 
#			## new code 
#			mdprime2 <- min(distance(X, Xorig)) 
#			if(mdprime2 < mdprime){ 
#				mdprime <- mdprime2 
#			}
#		} 
#		if(mdprime > md){ 	##accept
#			md <- mdprime	 
#		} else{ 
#			X[row,] <- xold ##reject
#		}
#	}
#	return(X) 
#}
##lhsRead = function()


##
#place = './modsSchnuteHHardFlatT30N150WWideExtra/'
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
#TT=30
#FtFmsy = rep(1, TT) #make faux catch
#
##
#B0 = 10000
#M  = 0.2
#
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
#???END
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


