rm(list=ls())

#
library(VGAM)
library(akima)
library(pracma)
library(mvtnorm)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)


#
#FUNCTIONS
#

#
streamlines <- function(x, y, u, v, step.dist=NULL, 
                        max.dist=NULL, col.ramp=c("white","black"), 
                        fade.col=NULL, length=0.05, ...) {

  ## Function for adding smoothed vector lines to a plot. 
  ## Interpolation powered by akima package

  ## step.distance - distance between interpolated locations (user coords)
  ## max.dist - maximum length of interpolated line (user coords)
  ## col.ramp - colours to be passed to colorRampPalette
  ## fade.col - NULL or colour to add fade effect to interpolated line
  ## ... - further arguments to pass to arrows

  ## build smoothed lines using interp function
  maxiter <- max.dist/step.dist
  l <- replicate(5, matrix(NA, length(x), maxiter), simplify=FALSE)
  names(l) <- c("x","y","u","v","col")
  l$x[,1] <- x
  l$y[,1] <- y
  l$u[,1] <- u
  l$v[,1] <- v
  for(i in seq(maxiter)[-1]) {
    l$x[,i] <- l$x[,i-1]+(l$u[,i-1]*step.dist)
    l$y[,i] <- l$y[,i-1]+(l$v[,i-1]*step.dist)
    r <- which(l$x[,i]==l$x[,i-1] & l$y[,i]==l$y[,i-1])
    l$x[r,i] <- NA
    l$y[r,i] <- NA
    for(j in seq(length(x))) {
      if(!is.na(l$x[j,i])) {
        l$u[j,i] <- c(interp(x, y, u, xo=l$x[j,i], yo=l$y[j,i])$z)
        l$v[j,i] <- c(interp(x, y, v, xo=l$x[j,i], yo=l$y[j,i])$z) 
      } 
    }
  }

  ## make colour a function of speed and fade line
  spd <- sqrt(l$u^2 + l$v^2) # speed
  spd <- apply(spd, 1, mean, na.rm=TRUE) # mean speed for each line
  spd.int <- seq(min(spd, na.rm=TRUE), max(spd, na.rm=TRUE), length.out=maxiter)
  cr <- colorRampPalette(col.ramp)
  cols <- as.numeric(cut(spd, spd.int))
  ncols <- max(cols, na.rm=TRUE)
  cols <- cr(ncols)[cols]
  if(is.null(fade.col)) {
    l$col <- replicate(maxiter, cols)
  } else {
    nfade <- apply(!is.na(l$x), 1, sum)
    for(j in seq(length(x))) {
      l$col[j,seq(nfade[j])] <- colorRampPalette(c(fade.col, cols[j]))(nfade[j])
    } 
  }

  ## draw arrows
  for(j in seq(length(x))) {
    arrows(l$x[j,], l$y[j,], c(l$x[j,-1], NA), c(l$y[j,-1], NA), 
           col=l$col[j,], length=0, ...)
    i <- which.max(which(!is.na(l$x[j,]))) # draw arrow at end of line
    if(i>1) {
      arrows(l$x[j,i-1], l$y[j,i-1], l$x[j,i], l$y[j,i], 
             col=l$col[j,i-1], length=length, ...) 
    }
  }

}

#
dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/(x+2))^2) }

#
getData = function(dir, xiSims, zetaSims){
        #dir    : a directory containing data
        #xiSims : xis of simulated data
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
        #X      :a vector of predictors
        #l      :a scalar length scale

        #
        sapply(X, function(x){K(x,X,l)})
}
#
Kn = function(X, Y, lx, ly){
        #X      :a vector of predictors
        #Y      :a vector of predictors
        #lx     :a scalar length scale for X
        #ly     :a scalar length scale for Y

        #
        sapply(Y, function(y){K(y,Y,ly)})*sapply(X, function(x){K(x,X,lx)})
}
#
Ku = function(X, Y, lx, ly){
        #X      :a vector of predictors
        #Y      :a vector of predictors
        #lx     :a scalar length scale for X
        #ly     :a scalar length scale for Y

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
zetaSims = rev(seq(0.1, 0.8, 0.01)) #rev(seq(0.25, 0.75, 0.01))
xiSims =   seq(0.5, 3.5, 0.05)

#
M = 0.2
#time: 70
D = getData(dir, xiSims, zetaSims)
D = D[D$cllhsV!=0,]
#D$cllhsV = 0

#
#GP
#

#
xiPred   = seq(0.3,3.7,length.out=30)      #seq(min(D$xi),max(D$xi),length.out=100)        #
zetaPred = seq(0.1,0.9,length.out=30)    #seq(min(D$zeta),max(D$zeta),length.out=100)    #
predMesh = expand.grid(xiPred, zetaPred)
colnames(predMesh) = c('xi', 'zeta')
#predMesh = meshgrid(xiPred, zetaPred)
#predMesh$X = as.vector(predMesh$X)
#predMesh$Y = as.vector(predMesh$Y)
##
#predMeshDF = as.data.frame(predMesh)
#colnames(predMeshDF) = c('xi', 'zeta')
line = lm(cllhs~xi+zeta, data=D)
pp = predict( line )
ppPred = predict(line, predMesh)
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
KPX = KnPred(predMesh$xi, predMesh$zeta, D$xi, D$zeta, lx, ly)
KXP = KnPred(D$xi, D$zeta, predMesh$xi, predMesh$zeta, lx, ly)
#KPP = KnPred(predMesh$xi, predMesh$zeta, predMesh$xi, predMesh$zeta, lx, ly)
#time: 66
ys = ppPred + KPX%*%KInv%*%(D$cllhs-pp) 

#
#OUTPUT
#

#
m = 3.5
bigD = mapply(function(xiHat, xi, zeta){
                dist(xiHat, xi, zeta)
        }, exp(ys)/M, xiPred, zetaPred
)
#
w = predMesh$xi>0.5 & predMesh$xi<3.5 & predMesh$zeta>0.2 & predMesh$zeta<0.75
u = (exp(ys)/M) - predMesh$xi  #exp(ys)/M #
u[!w] = 0 #u[u>(2)] = 0 #u[u<(-2)] = 0
v = (1/(exp(ys)/M+2)) - predMesh$zeta  #1/(exp(ys)/M+2) #
v[!w] = 0
#
pdf('gpBiasArrow.pdf')
#m = 3.5 #quantile(abs(bigD), 0.5)
#predMesh$xi, predMesh$zeta,
filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T),
        zlim=c(0, m),
        plot.axes = {
		quiver(predMesh$xi, predMesh$zeta, u, v)
               	#for(i in 1:length(predMesh$xi)){
		#	arrows(predMesh$xi[i], predMesh$zeta[i], predMesh$xi[i] + scale * u, predMesh$zeta + scale * v
		#}
               	lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3) 
        },
        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
)
dev.off()


#quiver(predMesh$xi, predMesh$zeta, u, v)



##
#pdf('streamers.pdf')
#plot(zeta~xi, data=predMesh, type="n", xlab="", ylab="")
#streamlines(predMesh$xi, predMesh$zeta, u, v, 
#		step.dist = 0.005, #1000, 
#		max.dist  = 0.1, #50000, 
#		col.ramp  = c("white","black"), 
#            	fade.col  = "white", 
#		length    = 0, 
#		lwd       = 2.5
#)
#lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#dev.off()
#
##
#m = 3.5
#
##
#bigD = mapply(function(xiHat, xi, zeta){
#                dist(xiHat, xi, zeta)
#        }, exp(ys)/M, xiPred, zetaPred
#)
##
#pdf('gpBiasStream.pdf')
##m = 3.5 #quantile(abs(bigD), 0.5)
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T),
#        zlim=c(0, m),
#        plot.axes = {
#		streamlines(predMesh$xi, predMesh$zeta, u, v, 
#			step.dist = 0.005, #1000, 
#			max.dist  = 0.1, #50000, 
#			col.ramp  = c("white","black"), 
#            		fade.col  = "white", 
#			length    = 0, 
#			lwd       = 2.5
#		)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#                #points(D$xi, D$zeta, pch='.')
#                #lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#                #points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#                #points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], p
#        },
#        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
#)
#dev.off()



#df$z[df$z<(-1.2)] = 0
#r = rasterFromXYZ(df, crs=CRS('+proj=longlat +datum=WGS84'))
#vectorplot(r, par.settings=RdBuTheme())




