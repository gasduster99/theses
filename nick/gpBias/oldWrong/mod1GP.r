rm(list=ls())

#
library(VGAM)
library(pracma)
library(mvtnorm)
library(plot.matrix)
library(RColorBrewer)
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
loglikeL = function(y, D, lx, ly){
	Kay = Kn(D$xi,D$zeta,lx,ly)+diag(D$cllhsV)
	return( dmvnorm(y,sigma=Kay, log=T) )
	
}

#
#DATA
#

#
dir = './modsFine/' #'./modsHess/'
##
#zetaSims = rev(c(0.15, seq(0.25, 0.75, 0.1)))
#xiSims = c(0.05, seq(0.5, 3.5, 0.5))
zetaSims = rev(seq(0.1, 0.80, 0.01)) #rev(seq(0.1, 0.8, 0.01))
xiSims =   seq(0.5, 3.5, 0.05)

#
M = 0.2
#time: 70
D = getData(dir, xiSims, zetaSims)
D = D[D$cllhsV!=0,]
#D$cllhsV = 0

#
xiPred   = seq(0.3,3.7,length.out=50)	#seq(0,4,length.out=100)	#seq(min(D$xi),max(D$xi),length.out=100)	#
zetaPred = seq(0.1,0.9,length.out=50)	#seq(0,0.8,length.out=100)	#seq(min(D$zeta),max(D$zeta),length.out=100)	#
predMesh = meshgrid(xiPred, zetaPred)
predMesh$X = as.vector(predMesh$X)
predMesh$Y = as.vector(predMesh$Y)
#
predMeshDF = as.data.frame(predMesh)
colnames(predMeshDF) = c('xi', 'zeta')
line = lm(cllhs~xi+zeta, data=D)
pp = predict( line )
ppPred = predict(line, predMeshDF)
#c(0.01, 0.001)
#c(0.0248, 0.00519)
#L = c(0.486, 0.121)
L = c(0.54, 0.134)
names(L) = c('lx', 'ly') 
#time: 340
t = system.time({ optOut = optim(L, function(x){-loglikeL(D$cllhs-pp, D, x[1], x[2])}) })
writeLines(sprintf("Opt: %s",t[3]))
lx = optOut$par['lx']  #sd(D$xi)/3
ly = optOut$par['ly']  #sd(D$zeta)/3
#
KInv = chol2inv(chol(Kn(D$xi,D$zeta,lx,ly)+diag(D$cllhsV)))
KPX = KnPred(predMesh$X, predMesh$Y, D$xi, D$zeta, lx, ly)
KXP = KnPred(D$xi, D$zeta, predMesh$X, predMesh$Y, lx, ly)
#KPP = KnPred(predMesh$X, predMesh$Y, predMesh$X, predMesh$Y, lx, ly)
#time: 66
ys = ppPred + KPX%*%KInv%*%(D$cllhs-pp) #cloglog(KPX%*%KInv%*%(D$cllhs-pp), inverse=T) #
##time: 470
#Ss = KPP-KPX%*%KInv%*%KXP
#writeLines(sprintf("Krig Var: %s",t[3]))

#
#CLLHS OUTPUT
#

#
pdf('gpCllhsFine.pdf')
filled.contour(xiPred, zetaPred, matrix(ys, nrow=length(xiPred), byrow=T), 
	zlim=c(-4, 4),#c(-3,0),
	plot.axes = {
		points(D$xi, D$zeta)
		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
		points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
	}
)
dev.off()

#
m = 3.5

#
pdf('gpXiBiasFine.pdf')
xDist = matrix((exp(ys)/M)-xiPred, nrow=length(xiPred), byrow=T)
#xDist[xDist<(-1.2)] = -10
filled.contour(xiPred, zetaPred, xDist, 
	zlim=c(-m, m),
	plot.axes = {
		points(D$xi, D$zeta, pch='.')
		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
	},
	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
)
dev.off()

#
pdf('gpZetaBiasFine.pdf')
yDist = matrix((1/(exp(ys)/M+2))-zetaPred, nrow=length(xiPred), byrow=T)
filled.contour(xiPred, zetaPred, yDist, 
	zlim=c(-0.55,0.55),
	plot.axes = {
		points(D$xi, D$zeta)
		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
	},
	color.palette = function(n) hcl.colors(n, "Blues")
)
dev.off()

#
bigD = mapply(function(xiHat, xi, zeta){ 
		dist(xiHat, xi, zeta) 
	}, exp(ys)/M, xiPred, zetaPred
) 
#
pdf('gpBiasFine.pdf')
#m = 3.5 #quantile(abs(bigD), 0.5)
filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T), 
	zlim=c(0, m),
	plot.axes = {
		points(D$xi, D$zeta, pch='.')
		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
	},
	color.palette = function(n) hcl.colors(n, "Blues", rev = TRUE)
)
dev.off()

#
w = predMesh$X>0.5 & predMesh$X<3.5 & predMesh$Y>0.2 & predMesh$Y<0.75
#
pdf('gpBiasArrow.pdf')
filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T),
        zlim=c(0, m),
        plot.axes = {
                quiver(predMesh$X[w[c(T,F)]], predMesh$Y[w[c(T,F)]], xDist[w[c(T,F)]], yDist[w[c(T,F)]])
                #for(i in 1:length(predMesh$xi)){
                #       arrows(predMesh$xi[i], predMesh$zeta[i], predMesh$xi[i] + scale * u, predMesh$zeta 
                #}
                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
        },
        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
)
dev.off()





#sapply(1:length(xiPred),
#	function(xi){sapply(1:length(zetaPred), 
#		function(zeta){ dist((exp(ys)/M)[], xiPred[], zetaPred[]) }
#	)
#}

##
#pdf('gpCllhsHessVar.pdf')
#filled.contour(xiPred, zetaPred, matrix(diag(Ss), nrow=length(xiPred), byrow=T), 
#	zlim=c(0, 0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#	}
#)
#dev.off()




##
#p = scatterplot3d(D$xi, D$zeta, D$cllhs)
#plane = lm(cllhs~xi+zeta, data=D)
#p$plane3d(plane)
#pp = predict(plane)
#scatterplot3d(D$xi, D$zeta, D$cllhs-pp)

##
##DISTANCE OUTPUT
##
#
##
#plane = lm(log(distHat-distMin)~xi+zeta, data=D)
#pp = predict(plane)
#yy = log(D$distHat-D$distMin)-pp
##
#yys = KPX%*%KInv%*%(yy)
##
#pdf('gpLogBiasNoHess.pdf')
#filled.contour(xiPred, zetaPred, matrix(yys, nrow=length(xiPred), byrow=T),
#        #zlim=c(-50, 20),
#        plot.axes = {
#                points(D$xi, D$zeta)
#                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3) 
#                #points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#                #points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#        }
#)
#dev.off()
