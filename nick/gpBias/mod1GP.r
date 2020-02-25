rm(list=ls())

#
library(VGAM)
library(pracma)
library(plot.matrix)
library(scatterplot3d)


#
#FUNCTIONS
#

#
dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/(x+2))^2) }

#
getData = function(dir, xiSims, zetaSims){
	#dir	: a directory containing data
	#xiSims	: xis of simulated data
	#zetaSims: zetas of simulated data
	
	#
	di = 1
	D = data.frame(xi=double(), zeta=double(), xiHat=double(), xiMin=double(), zetaHat=double(), distHat=double(), distMin=double(), cllhs=double(), cllhsV=double(), stringsAsFactors=F)
	for(i in 1:length(zetaSims)){
        	for(j in 1:length(xiSims)){
                	#
                	fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
                	fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
			
			#
			if( file.exists(fileDat) & file.exists(fileFit) ){
                		#
                		dat = readRDS(fileDat)
                		fit = readRDS(fileFit)

                		#
                		xiHat   = -log(1-cloglog(fit$cllhs, inverse=T))/M
                		zetaHat = 1/(xiHat+2)
                		distHat = sqrt((dat$xi-xiHat)^2 + (dat$zeta-zetaHat)^2)
                		opt = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)
                		#
				if(length(fit$rsCov)==0){ v=0 }else{ v=fit$rsCov['cllhs', 'cllhs'] }
                		D[di,]  = c(dat$xi, dat$zeta, xiHat, opt$minimum, zetaHat, distHat, opt$objective, fit$cllhs, v)
                		di=di+1
                	}
		}
	}
	#
	return(D)
}

#all parameters scalar
K = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
Kmat = function(X, l){ 
	#X	:a vector of predictors
	#l	:a scalar length scale
	
	#
	sapply(X, function(x){K(x,X,l)})
}
#
Kn = function(X, Y, lx, ly){ 
	#X	:a vector of predictors
	#Y	:a vector of predictors
	#lx	:a scalar length scale for X
	#ly	:a scalar length scale for Y
	
	#
	sapply(Y, function(y){K(y,Y,ly)})*sapply(X, function(x){K(x,X,lx)})
}
#
Ku = function(X, Y, lx, ly){ 
	#X	:a vector of predictors
	#Y	:a vector of predictors
	#lx	:a scalar length scale for X
	#ly	:a scalar length scale for Y
	
	#
	sapply(Y, function(y){K(y,Y,ly)})+sapply(X, function(x){K(x,X,lx)})
}
#
KnPred = function(Xs, Ys, X, Y, lx, ly){
	sapply(Y, function(y){K(y,Ys,ly)})*sapply(X, function(x){K(x,Xs,lx)})
}

#
#DATA
#

#
dir = './modsHess/'
#
zetaSims = rev(seq(0.25, 0.75, 0.1))
xiSims = seq(0.5, 3.5, 0.5)

#
M = 0.2
#
D = getData(dir, xiSims, zetaSims)
#D$cllhsV = 0

#
xiPred   = seq(0,4,length.out=100)	#seq(min(D$xi),max(D$xi),length.out=100)	#
zetaPred = seq(0,1,length.out=100)	#seq(min(D$zeta),max(D$zeta),length.out=100)	#
predMesh = meshgrid(xiPred, zetaPred)
predMesh$X = as.vector(predMesh$X)
predMesh$Y = as.vector(predMesh$Y)
#
lx = 1/2*sd(D$xi)
ly = 1/2*sd(D$zeta)
KInv = chol2inv(chol(Kn(D$xi,D$zeta,lx,ly)+diag(D$cllhsV)))
KPX = KnPred(predMesh$X, predMesh$Y, D$xi, D$zeta, lx, ly)
KXP = KnPred(D$xi, D$zeta, predMesh$X, predMesh$Y, lx, ly)
KPP = KnPred(predMesh$X, predMesh$Y, predMesh$X, predMesh$Y, lx, ly)
#
pp = predict( lm(cllhs~xi+zeta, data=D) )
ys = KPX%*%KInv%*%(D$cllhs-pp) #cloglog(KPX%*%KInv%*%(D$cllhs-pp), inverse=T) #
#Ss = KPP-KPX%*%KInv%*%KXP

#
#OUTPUT
#

#
filled.contour(xiPred, zetaPred, matrix(cloglog(ys, inverse=T), nrow=length(xiPred), byrow=T), 
	plot.axes = {
		points(D$xi, D$zeta)
		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
		points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
		points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
	}
)

##
#dev.new()
#filled.contour(xiPred, zetaPred, matrix(diag(Ss), nrow=length(xiPred), byrow=T), 
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#	}
#)

##
#p = scatterplot3d(D$xi, D$zeta, D$cllhs)
#plane = lm(cllhs~xi+zeta, data=D)
#p$plane3d(plane)
#pp = predict(plane)
#scatterplot3d(D$xi, D$zeta, D$cllhs-pp)



