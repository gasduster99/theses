rm(list=ls())

#
library(VGAM)
library(boot)
library(pracma)
library(mvtnorm)
library(graphics)
library(parallel)
library(rootSolve)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)
#
source("gpClass0.0.1.r")

#
#FUNCTIONS
#

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
#P0Funk = function(alpha, beta, gamma, M){ (alpha/(M)-1)^gamma * beta^-gamma }
#exp(gamma*log(alpha/(M+ff)-1) - gamma*log(beta)) }
PBar = function(ff, alpha, beta, gamma, M){ (alpha/(M+ff)-1)^gamma * beta^-gamma }

#
getBeta = function(alpha, gamma, M, P0){
        optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-P0)^2 }, method='Brent', lower=eps(), upper=10^6)$par
}

#
FMsy = function(alpha, gamma, M){
        #
        FUpper = alpha-M
        root = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
        #
        if(length(root)<1){ return(NA) }
        if(length(root)>1){ return(root[1]) }
        #
        return(root)
}

#
dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/(x+2))^2) }

##
#getData = function(dir, xiSims, zetaSims){
#	#dir	: a directory containing data
#	#xiSims	: xis of simulated data
#	#zetaSims: zetas of simulated data
#	
#	#
#	di = 1
#	D = data.frame(xi=double(), zeta=double(), xiHat=double(), xiMin=double(), zetaHat=double(), distHat=double(), distMin=double(), cllhs=double(), cllhsV=double(), stringsAsFactors=F)
#	for(i in 1:length(zetaSims)){
#        	for(j in 1:length(xiSims)){
#                	#
#                	fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
#                	fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
#			
#			#
#			if( file.exists(fileDat) & file.exists(fileFit) ){
#                		#
#                		dat = readRDS(fileDat)
#                		fit = readRDS(fileFit)
#
#               		#
#                		xiHat   = -log(1-cloglog(fit$cllhs, inverse=T))/M
#                		zetaHat = 1/(xiHat+2)
#                		distHat = sqrt((dat$xi-xiHat)^2 + (dat$zeta-zetaHat)^2)
#                		opt = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)
#                		#
#				if(length(fit$rsCov)==0){ v=0 }else{ v=fit$rsCov['cllhs', 'cllhs'] }
#                		D[di,]  = c(dat$xi, dat$zeta, xiHat, opt$minimum, zetaHat, distHat, opt$objective, fit$cllhs, v)
#                		di=di+1
#                	}
#		}
#	}
#	#
#	return(D)
#}

#
getData = function(dir, xiSims, zetaSims){
	#dir	: a directory containing data
	#xiSims	: xis of simulated data
	#zetaSims: zetas of simulated data
	
	#
	di = 1
	D = data.frame(xi=double(), zeta=double(), xiHat=double(), zetaHat=double(), s2=double(), minDist=double(), stringsAsFactors=F)
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
				xiHat   = FMsy(fit$alpha, fit$gamma, M)/M 
				zetaHat = PBar(xiHat*M, fit$alpha, fit$beta, fit$gamma, M)/PBar(0, fit$alpha, fit$beta, fit$gamma, M)
				md = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)$objective				

				#NOTE: replace with a propper observation uncertainty
				if(length(fit$rsCov)==0){ v=0 }else{ v=fit$rsCov['lalpha', 'lalpha'] }
				D[di,]  = c(dat$xi, dat$zeta, xiHat, zetaHat, v, md)
				di=di+1
			}
		}
	}
	#
	return(D)
}

#X's are 2-vectors; L is 2x2; s2 is scalar
S2 = function(X0, X1, s2, v0, v1, cr){
	maxD = dbinorm(X0[,1], X0[,2], X0[,1], X0[,2], v0, v1, sqrt(v0*v1)*cr, log=T)
	s2*mcmapply(function(x01, x02, m){
		exp(dbinorm(X1[,1], X1[,2], x01, x02, v0, v1, sqrt(v0*v1)*cr, log=T)-m)
	}, X0[,1], X0[,2], maxD, mc.cores=detectCores())	
}

#
logitM = function(n, mu, sig, K=10^5){
	mean( inv.logit(qnorm((1:(K-1))/K, mu, sig))^n )
}

#
#DATA
#

#
dir = "./modsShepTry/"#'./modsFine/' #'./modsHess/'
#
zetaSims = seq(0.1, 0.8, 0.05) 	#rev(seq(0.1, 0.80, 0.01)) #
xiSims =   seq(0.5, 3.5, 0.25)		#seq(0.5, 3.5, 0.05)       #

#
M = 0.2
#time: 70
D = getData(dir, xiSims, zetaSims)
D = D[D$s2!=0 & D$xiHat<20,]
#
png('shepDat.png')
plot(D[1:2], 
	ylim=c(min(D[c(2,4)]), max(D[c(2,4)])), 
	xlim=c(min(D[c(1,3)]), max(D[c(1,3)])),
	main='Shephard',
	col=map2color(D$minDist, hcl.colors(60, "Zissou 1", rev=T))
)
points(D[3:4], col=map2color(D$minDist, hcl.colors(60, "Zissou 1", rev=T)))
dev.off()

#X = cbind(D$xi, D$zeta)
#
##
#xiPred   = seq(0.7,3.5,length.out=50)	#seq(0.5, 3.5, 0.01) 		#	
#zetaPred = seq(0.2,0.75,length.out=50)	#seq(0.1, 0.75, 0.01)   	#
#predMesh = as.matrix(expand.grid(xiPred, zetaPred))
#colnames(predMesh) = c('xi', 'zeta')
##
#line = lm(cllhs~xi+zeta, data=D)
##gp = gpModel$new( S2=S2, X=X, Y=D$cllhs, obsV=D$cllhsV, #lm=line, 
##        B  = line$coeff,	#c(4.1894, -0.1452, -10.3066), 
##	v0 = 0.343979394716728,
##        v1 = 0.0270240223433852,
##        s2 = 2.93501281,
##        cr = -0.433574770896697
##)
###
##optOut0 = gp$fit(c('v0', 'v1', 's2', 'cr'),
##        lower   = c(eps(), eps(), eps(), -1),
##        upper   = c(Inf, Inf, Inf, 1),
##        cov     = T
##)
##gp$save('gpHotStart.rda')
#gp = readRDS('gpHotStart.rda')
#gp$printSelf()
##
#bnd = rep(Inf, length(line$coeff))
#optOut1 = gp$fit(c('B'),
#        lower   = c(-bnd),
#        upper   = c(bnd),
#        cov     = T
#)
#gp$save('gpHotterStart.rda')
#gp$printSelf()
##
#bnd = rep(Inf, length(line$coeff))
#optOut2 = gp$fit(c('v0', 'v1', 's2', 'cr', 'B'),
#        lower   = c(eps(), eps(), eps(), -1, -bnd),
#        upper   = c(Inf, Inf, Inf, 1, bnd),
#	#gaBoost = list(run=30, popSize=detectCores(), parallel=detectCores()),
#        cov     = T
#)
#gp$save('gpBigOpt.rda')
#gp$printSelf()
#
##
#yy = gp$predictMean(predMesh)
##SS = diag(gp$predictVar(predMesh))
###
##xS = exp(2*yy+SS)*(exp(SS)-1)/M^2
##K = 10^6
##yS = mcmapply(function(yy, SS){
##		logitM(2, log(2*M)-yy, SS, K=K) - logitM(1, log(2*M)-yy, SS, K=K)^2
##	}, yy, SS, mc.cores=8
##)
##xS = matrix(xS, nrow=length(xiPred))
##yS = matrix(yS, nrow=length(xiPred))
#
##
#xDist = matrix((exp(yy)/M)-xiPred, nrow=length(xiPred)) #, byrow=T)
#yDist = matrix((1/(exp(yy)/M+2))-zetaPred, nrow=length(xiPred))#, byrow=T)
##
#bigD = mcmapply(function(xiHat, xi, zeta){
#                dist(xiHat, xi, zeta)
#        }, exp(yy)/M, xiPred, zetaPred, mc.cores=detectCores()
#)
#
##
##OUTPUT
##
#
##
#m = 2.5
#
##
#pdf('gpXiBiasFineRB.pdf')
##dev.new()
#filled.contour(xiPred, zetaPred, xDist, 
#	zlim=c(-m, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#		#
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#	},
#	plot.title = title(
#		main = "Bias in Estimated Optimal Fishing Rate",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#pdf('gpZetaBiasFineRB.pdf')
##dev.new()
#filled.contour(xiPred, zetaPred, yDist, 
#	zlim=c(-0.55,0.55),
#        ylab = expression(frac(F[msy]/F[0])),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')	
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#	},
#	plot.title = title(
#		main = "Bias in Estimated Optimal Biomass",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3"),
#)
#dev.off()
#
##
#pdf('gpBiasArrowB.pdf')
##dev.new()
#datGrid = expand.grid(xiSims, zetaSims)
#colnames(datGrid) = c('xi', 'zeta')
#w = (predMesh[,'xi']>0.5 & predMesh[,'xi']<3.5 & predMesh[,'zeta']>0.2 & predMesh[,'zeta']<0.75) #paste0(predMesh[,'xi'], predMesh[,'zeta'])%in%paste0(datGrid[,'xi'], datGrid[,'zeta']) #
#thin = c(T,rep(F,3))
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred)),
#        zlim=c(0, m),
#        plot.axes = {
#                quiver(
#			predMesh[w,'xi'][thin], predMesh[w,'zeta'][thin],
#			xDist[w][thin], yDist[w][thin],
#			scale=0.05
#		) 
#                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#        },
#	plot.title = title(
#		main = "Directional Bias",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	#key.title = title(main = "Directional Bias"),
#        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
#)
#dev.off()
















#
#dev.new()
#filled.contour(xiPred, zetaPred, yDist/sqrt(yS), 
#	zlim=c(-0.6,0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)




##all parameters scalar
#Kr = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
#K = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
#Kmat = function(X, l){ 
#	#X	:a vector of predictors
#	#l	:a scalar length scale
#	
#	#
#	sapply(X, function(x){K(x,X,l)})
#}
##
#Kn = function(X, Y, lx, ly){ 
#	#X	:a vector of predictors
#	#Y	:a vector of predictors
#	#lx	:a scalar length scale for X
#	#ly	:a scalar length scale for Y
#	
#	#
#	sapply(Y, function(y){K(y,Y,ly)})*sapply(X, function(x){K(x,X,lx)})
#}
##
#Ku = function(X, Y, lx, ly){ 
#	#X	:a vector of predictors
#	#Y	:a vector of predictors
#	#lx	:a scalar length scale for X
#	#ly	:a scalar length scale for Y
#	
#	#
#	sapply(Y, function(y){K(y,Y,ly)})+sapply(X, function(x){K(x,X,lx)})
#}
##
#KnPred = function(Xs, Ys, X, Y, lx, ly){
#	sapply(Y, function(y){K(y,Ys,ly)})*sapply(X, function(x){K(x,Xs,lx)})
#}
#
##
#loglikeL = function(y, D, lx, ly, r, s2){	#function(y, D, L, sig2){ #
#	#
#	X = cbind(D$xi, D$zeta)
#	SIG2 = S2(X, X, lx, ly, r, s2) + diag(D$cllhsV)
#		
#	#
#	return( dmvnorm(y, sigma=SIG2, log=T) )
#}





#pp = predict( line )
#ppPred = predict(line, predMesh)


###
##optOut = c(mean(c(0.34395317, 0.02705277)), 2.93501281) 
##names(optOut) = c('lx', 'sig2')
#optOut = c(0.34395317, 0.02705277, -0.43195056, 2.93501281) 
#names(optOut) = c('lx', 'ly', 'cr', 'sig2')
##optOut = c(0.34395317, 0.02705277, 2.93501281)
##names(optOut) = c('lx', 'ly', 'sig2') 
##time: 340
#t = system.time({ 
#	optOut = optim(optOut, 
#		function(x){ 	
#			-loglikeL(D$cllhs-pp, D, x[1], x[2], x[3], x[4])
#		},
#		lower = c(10^-6, 10^-6, -1, 0),
#		upper = c(Inf, Inf, 1, Inf),
#		method = "L-BFGS-B",
#		hessian = T
#	)
#})
#writeLines(sprintf("Opt: %s",t[3]))
##IL = matrix(c(optOut$par[1], optOut$par[3], optOut$par[3], optOut$par[2]), nr=2)
#lx   = optOut$par['lx']
#ly   = optOut$par['ly']
#r    = optOut$par['cr']
#sig2 = optOut$par['sig2']
## 
##KInv = chol2inv(chol( sig2*Kn(D$xi,D$zeta,lx,ly) + diag(D$cllhsV) ))
#X = cbind(D$xi, D$zeta)
#KInv = chol2inv(chol( S2(X, X, lx, ly, r, sig2) + diag(D$cllhsV) ))
##KPX  = sig2*KnPred(predMesh$xi, predMesh$zeta, D$xi, D$zeta, lx, ly)
#KPX = S2(X, predMesh, lx, ly, r, sig2)
##KXP  = sig2*KnPred(D$xi, D$zeta, predMesh$xi, predMesh$zeta, lx, ly)
##KPP  = sig2*KnPred(predMesh$xi, predMesh$zeta, predMesh$xi, predMesh$zeta, lx, ly)
####time: 66
#ys = ppPred + KPX%*%KInv%*%(D$cllhs-pp) #cloglog(KPX%*%KInv%*%(D$cllhs-pp), inverse=T) #
###time: 470
##Ss = KPP-KPX%*%KInv%*%KXP
###writeLines(sprintf("Krig Var: %s",t[3]))

##
##CLLHS OUTPUT
##
#
##
#pdf('gpCllhsFineR.pdf')
#filled.contour(xiPred, zetaPred, matrix(ys, nrow=length(xiPred), byrow=T), 
#	zlim=c(-3.5, 1.5), #c(-3,0),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	}
#)
#dev.off()
#
##
#m = 3.5
#
##
#pdf('gpXiBiasFineR.pdf')
#xDist = matrix((exp(ys)/M)-xiPred, nrow=length(xiPred), byrow=T)
##xDist[xDist<(-1.2)] = -10
#filled.contour(xiPred, zetaPred, xDist, 
#	zlim=c(-m, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#pdf('gpZetaBiasFineR.pdf')
#yDist = matrix((1/(exp(ys)/M+2))-zetaPred, nrow=length(xiPred), byrow=T)
#filled.contour(xiPred, zetaPred, yDist, 
#	zlim=c(-0.6,0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#bigD = mapply(function(xiHat, xi, zeta){ 
#		dist(xiHat, xi, zeta) 
#	}, exp(ys)/M, xiPred, zetaPred
#) 
##
#pdf('gpBiasFineR.pdf')
##m = 3.5 #quantile(abs(bigD), 0.5)
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T), 
#	zlim=c(0, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blues", rev = TRUE)
#)
#dev.off()


###
##w = predMesh$X>0.5 & predMesh$X<3.5 & predMesh$Y>0.2 & predMesh$Y<0.75
###
##pdf('gpBiasArrow.pdf')
##filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T),
##        zlim=c(0, m),
##        plot.axes = {
##                quiver(predMesh$X[w[c(T,F)]], predMesh$Y[w[c(T,F)]], xDist[w[c(T,F)]], yDist[w[c(T,F)]])
##                #for(i in 1:length(predMesh$xi)){
##                #       arrows(predMesh$xi[i], predMesh$zeta[i], predMesh$xi[i] + scale * u, predMesh$zeta 
##                #}
##                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
##        },
##        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
##)
##dev.off()
#
#
#
#
#
##sapply(1:length(xiPred),
##	function(xi){sapply(1:length(zetaPred), 
##		function(zeta){ dist((exp(ys)/M)[], xiPred[], zetaPred[]) }
##	)
##}
#
###
##pdf('gpCllhsHessVar.pdf')
##filled.contour(xiPred, zetaPred, matrix(diag(Ss), nrow=length(xiPred), byrow=T), 
##	zlim=c(0, 0.6),
##	plot.axes = {
##		points(D$xi, D$zeta)
##		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
##	}
##)
##dev.off()
#
#
#
#
###
##p = scatterplot3d(D$xi, D$zeta, D$cllhs)
##plane = lm(cllhs~xi+zeta, data=D)
##p$plane3d(plane)
##pp = predict(plane)
##scatterplot3d(D$xi, D$zeta, D$cllhs-pp)
#
###
###DISTANCE OUTPUT
###
##
###
##plane = lm(log(distHat-distMin)~xi+zeta, data=D)
##pp = predict(plane)
##yy = log(D$distHat-D$distMin)-pp
###
##yys = KPX%*%KInv%*%(yy)
###
##pdf('gpLogBiasNoHess.pdf')
##filled.contour(xiPred, zetaPred, matrix(yys, nrow=length(xiPred), byrow=T),
##        #zlim=c(-50, 20),
##        plot.axes = {
##                points(D$xi, D$zeta)
##                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3) 
##                #points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
##                #points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
##        }
##)
##dev.off()
