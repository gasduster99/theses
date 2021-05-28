rm(list=ls())

#
library(plgp)
library(boot)
library(mlegp)
library(pracma)

#
#FUNCTIONS
#

#
maximin = function(n, m, TT=100000) {
	#n : number of points to draw
	#m : dimension of space to draw from
	#
	X = matrix(runif(n*m), ncol=m) 
	d = distance(X)
	d = d[upper.tri(d)]
	#
	md = min(d)
	for(t in 1:TT) {
		#propse a row
		row = sample(1:n, 1)
		xold = X[row,]
		#assume we accept
		X[row,] = runif(m)
		d = distance(X) #new distance
		d = d[upper.tri(d)]
		mdprime = min(d)
		if( mdprime > md ){ md=mdprime ## accept 
		} else { X[row,] <- xold } ## reject
	}
	return(X) 
}

#
maximinAdd = function(n, m, TT=100000, XStart=NULL) {
	#n : number of points to draw
	#m : dimension of space to draw from
	#
	#value :  only retunrs new additions
	#
	X = matrix(runif(n*m), ncol=m) 
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
		X[row,] = runif(m)
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
maximinF = function(n, m, f, Xrange, Yrange, TT=100000, XStart=NULL) {
	#n : number of points to draw
	#m : dimension of space to draw from
	#f : function to spread over
	#Xrange : limits in X
	# I need to figure out how to actually place the limits on the space of y=f(x)'s
	#Yrange : limits in f(X) 
	#
	#value :  only retunrs new additions
	
	#
	Xrange = matrix(Xrange, nrow=m)
	X = matrix(rep(NA, n*m), ncol=m)
	#
	jobsNotDone = rep(T, n)
	while( any(jobsNotDone) ){
		#
		for(i in (1:n)[jobsNotDone]){
			#
			X[i,]=runif(m, Xrange[,1], Xrange[,2])
			#for(j in 1:m){ 
			#	X[i,j]=runif(1, Xrange[j,1], Xrange[j,2]) 
			#}
		}
		FF = f(X)
		jobsNotDone = !(FF>Yrange[1] & FF<Yrange[2]) 
	}
	# 
	d = distance(FF)
	d = d[upper.tri(d)]
	#
	md = min(d)
	#handle initiated design
	if(!is.null(XStart)){
		FFStart = f(XStart)
		md2 = min(distance(FF, FFStart))
		if(md2 < md){ md=md2 }
	}
	for(t in 1:TT){
		#propose a row
		row = sample(1:n, 1)
		xold = X[row,]
		#
		jobNotDone = T 	
		while( jobNotDone ){
			#assume we accept 
			X[row,]=runif(m, Xrange[,1], Xrange[,2])
			#for(j in 1:m){ X[row,j]=runif(1, Xrange[j,1], Xrange[j,2]) }
			FF = f(X)	
			jobNotDone = !(FF[row]>Yrange[1] & FF[row]<Yrange[2]) 
		}
		d = distance(FF) #new distance
		d = d[upper.tri(d)]
		mdprime = min(d)
		#handle initial design (make min dist smaller if neccesary with DStart included)
		if(!is.null(XStart)){
			mdprime2 = min(distance(FF, FFStart))
			if(mdprime2 < mdprime) mdprime=mdprime2 
		}
		#
		if( mdprime > md ){ md=mdprime ## accept 
		} else{ X[row,]=xold } ## reject
	}
	#return only the augmenting points
	return(X) 
}

##
##GP
##
#
##
#M = 10^4
#starting = list("phi"=1, "sigma.sq"=0.25)#, "tau.sq"=1, "nu"=1)
#tuning   = list("phi"=0.1, "sigma.sq"=0.1)#, "tau.sq"=0.1, "nu"=1)
#priors   = list(
#        "beta.Flat",
#        "phi.Unif"    = c(0.3, 5),
#        "sigma.sq.IG" = c(2, 0.1)#,
#        #"tau.sq.IG"   = c(2, 0.01),
#        #"nu.unif"     = c(0.1, 10)
#)
#cov.model = "gaussian" #"matern" #"exponential"
#spB = spLM(y~X, 
#       coords    = cbind(1,X),
#       #knots     = c(8, 6, 0.25), #c(16, 12, 0.25), #c(7, 5, -0.1), #  #turn on the predictive process
#       n.samples = M,
#       cov.model = cov.model,
#       priors    = priors,  
#       starting  = starting,  
#       tuning    = tuning,
#       #n.report  = M*0.1,
#       #n.omp.threads=8,
#       verbose   = F
#)
##
#capture.output(
#spBPred <- spPredict(spB, cbind(1,x), cbind(1,x), 
#       start   = 200, 
#       thin    = 8, 
#       verbose = F 
#       #n.report= ((M-200)/8+1)*0.1, 
#       #n.omp.threads=8
#), file="/dev/null")
##
#return( rowMeans(spBPred$p.y.predictive.samples, na.rm=T) )


#
#MAIN
#

#
f = function(x) inv.logit(x) 
Xlimits = c(-10, 10)
Ylimits = c(0, 1)
##
#f = function(x) sin(x)+ x/2
#Xlimits = c(-6, 6) #c(-pi/2, pi/2)
#Ylimits = f(Xlimits)
##
XStar = seq(Xlimits[1], Xlimits[2], length.out=5*(Xlimits[2]-Xlimits[1]))

#f will be made of GP model prediction
nStart = 10
#X = maximinF(nStart, 1, f, Xlimits, Ylimits, TT=100000, XStart=NULL)
X = runif(nStart, Xlimits[1], Xlimits[2])
y = f(X)
#
Xall = t(t(X))
Yall = t(t(y))

#
curve(f, Xlimits[1], Xlimits[2], lwd=3)
points(X, f(X))
rug(X)
points(rep(Xlimits[1], length(X)), y, pch='-', cex=2)

#
M = 100-nStart
for(i in 1:M){
	#FIT and update stat model
	gpFit = mlegp(Xall, Yall, constantMean=F, verbose=0)
	fHat = function(x, d=1){
		#
		x = matrix(x, ncol=d)
		return( predict(gpFit, x) )
	}
	
	#plot new predicitons
	lines(XStar, fHat(XStar), col=i)	

	#find new points to evaluate
	X = maximinF(1, 1, fHat, Xlimits, Ylimits, TT=25000, XStart=Xall)
	y = f(X)
	#add nnew evaluations to all evaluations
	Xall = rbind(Xall, X)
	Yall = rbind(Yall, y)

	#plot new evaluations	
	text(X, y, label=as.character(nStart+i), col=i)
	text(X, rep(Ylimits[1], length(y)), label=as.character(nStart+i), col=i)
	text(rep(Xlimits[1], length(X)), y, label=as.character(nStart+i), col=i)
}


##f is updated with new GP model predictor and cycle
#X2 = maximinF(10, 1, fHat, Xlimits, Ylimits, TT=10000, XStart=X)
#XX2 = rbind(X, X2)
#yy2 = f(XX2)
#
##
#gpFit = mlegp(XX2, yy2, constantMean=F)
#
##
#points(X2, f(X2), col='blue')
#rug(X2, col='blue')
#points(rep(Xlimits[1], length(X2)), f(X2), pch='-', col='blue')
#lines(XStar, fHat(XStar), col='red', lwd=2)
#
##
#X3 = maximinF(10, 1, fHat, Xlimits, Ylimits, TT=100000, XStart=rbind(X, X2))
##
#XX3 = rbind(X, X2, X3)
#yy3 = f(XX3)
#gpFit = mlegp(XX3, yy3, constantMean=F)
#
##
#points(X3, f(X3), col='red')
#rug(X3, col='red')
#points(rep(Xlimits[1], length(X3)), f(X3), pch='-', col='red')
#lines(XStar, fHat(XStar), col='green')



