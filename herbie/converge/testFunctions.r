figPlace = "~/Documents/school/ucscGrad/thesis/herbie/writting/figures/"

colPersp = function(x,y,z,pal,nb.col,...,xlg=TRUE,ylg=TRUE){
	colnames(z) <- y
	rownames(z) <- x

	nrz <- nrow(z)
	ncz <- ncol(z) 

	color <- pal(nb.col)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	facetcol <- cut(zfacet, nb.col)
	par(xlog=xlg,ylog=ylg)
	persp(
		as.numeric(rownames(z)),
		as.numeric(colnames(z)),
		as.matrix(z),
		col=color[facetcol],
		...
		)
}

#x shall have grid spaces in the rows and dimensions in the columns
rosenbrock = function(x){
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

rastrigin = function(x){
        out = matrix(NA, nrow=dim(x)[1], ncol=1)
        for (i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = 10*length(ex)+sum(ex^2-10*cos(2*pi*ex))
        }
        return(out)
}

easom = function(x){
	out = -apply(cos(x), 1, prod) * exp(-apply((x-pi)^2, 1, sum))
	#out = matrix(NA, nrow=dim(x)[1], ncol=1)	
	#for(i in seq(1, dim(x)[1])){
	#	ex = x[i,]
	#	out[i] = -prod(cos(ex)) * exp( -sum((ex-pi)^2) )
	#}
	return(out)
}


#
#ROSENBROCK
#
h=0.01
dimX = seq(-2, 2, h)
dimY = seq(-3, 5, h)
#
Z = matrix(NA, nrow=length(dimX), ncol=length(dimY))
logZ = matrix(NA, nrow=length(dimX), ncol=length(dimY))
for(i in seq(1, length(dimX))){
	for(j in seq(1, length(dimY))){
		Z[i,j] = rosenbrock( matrix(c(dimX[i], dimY[j]), nrow=1, ncol=2) )#shekel(c(dimD[i], dimD[j]), m=9)#easom(matrix(c(dimD[i], dimD[j]), ncol=2)) 
		logZ[i,j] = log(Z[i,j])
	}
}
##
#pdf(sprintf("%sroseContour.pdf", figPlace))
## define jet colormap
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#filled.contour( dimX, dimY, Z,
#                        color.palette = jet.colors,
#                        xlab=expression(x[1]),
#                        ylab=expression(x[2]),
#                        #main='Rosenbrock',
#			#zlim=c(0, 3000)
#			plot.axes={ points(c(1), c(1), col="white", pch=20); axis(1); axis(2) }
#)
#dev.off()
#
pdf(sprintf("%srosePersp.pdf", figPlace))
colPersp( dimX, dimY, Z, 
	jet.colors, 
	200, 
	phi=45, 
	theta=-40, 
	border=NA,
	box=T,
	ticktype="detailed",
	xlab="x1",
	ylab="x2",
	#zlab='Rosenbrock'
	zlab=""
	#main='Rosenbrock'
)
dev.off()
#

##
##RASTRIGIN
##
#h=0.01
#dimX = seq(-2.5, 2.5, h)
#dimY = seq(-2.5, 2.5, h)
##
#Z = matrix(NA, nrow=length(dimX), ncol=length(dimY))
#logZ = matrix(NA, nrow=length(dimX), ncol=length(dimY))
#for(i in seq(1, length(dimX))){
#	for(j in seq(1, length(dimY))){
#		Z[i,j] = rastrigin( matrix(c(dimX[i], dimY[j]), nrow=1, ncol=2) )#shekel(c(dimD[i], dimD[j]), m=9)#easom(matrix(c(dimD[i], dimD[j]), ncol=2)) 
#		logZ[i,j] = log(Z[i,j])
#	}
#}
##
##dev.new()
#pdf(sprintf("%srastContour.pdf", figPlace))
## define jet colormap
#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#filled.contour( dimX, dimY, Z,
#                        color.palette = jet.colors,
#                        xlab=expression(x[1]),
#                        ylab=expression(x[2]),
#                        #main='Rastrigin',
#			#zlim=c(0, 3000)
#			plot.axes={ points(c(0), c(0), col="white", pch=20); axis(1); axis(2) }			
#)
#dev.off()
#
##dev.new()
#pdf(sprintf("%srastPersp.pdf", figPlace))
#colPersp( dimX, dimY, Z, 
#	jet.colors, 
#	200, 
#	phi=45, 
#	theta=-40, 
#	border=NA,
#	box=T,
#	ticktype="detailed",
#	xlab="x1",
#	ylab="x2",
#	#zlab='Rosenbrock'
#	zlab=""#,
#	#main='Rastrigin'
#)
#dev.off()
##

##
##EASOM
##
#h=0.025
#dimX = seq(-11, 11, h)
#dimY = seq(-11, 11, h)
##
#Z = matrix(NA, nrow=length(dimX), ncol=length(dimY))
#logZ = matrix(NA, nrow=length(dimX), ncol=length(dimY))
#for(i in seq(1, length(dimX))){
#	for(j in seq(1, length(dimY))){
#		Z[i,j] = easom( matrix(c(dimX[i], dimY[j]), nrow=1, ncol=2) )
#	}
#}
##
##dev.new()
#pdf(sprintf("%seasomContour.pdf", figPlace))
## define jet colormap
#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#222222"))
#filled.contour( dimX, dimY, Z,
#                        color.palette = jet.colors,
#                        xlab=expression(x[1]),
#                        ylab=expression(x[2]),
#                        #main='Easom',
#			#zlim=c(0, 3000)
#			plot.axes={ points(c(pi), c(pi), col="white", pch=20); axis(1); axis(2) }
#)
#dev.off()
#
##dev.new()
#pdf(sprintf("%seasomPersp.pdf", figPlace))
#colPersp( dimX, dimY, Z, 
#	jet.colors, 
#	200, 
#	phi=35, 
#	theta=180, 
#	border=NA,
#	box=T,
#	ticktype="detailed",
#	xlab="x1",
#	ylab="x2",
#	#zlab='Rosenbrock'
#	zlab=""#,
#	#main='Easom'
#)
#dev.off()
##
