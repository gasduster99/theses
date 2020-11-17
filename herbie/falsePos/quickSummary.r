rm(list=ls())

#
library(GGally)
library(pracma)
library(grDevices)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)

#
#FUNCTIONS
#

#
unpackIts = function(out){
	#
	outter = c()#data.frame(lamda=numeric(), W=integer(), it=integer())
	for(i in 1:length(out)){
		#
		its = unlist(out[[i]]$itConv[-(1:(length(wGrid)+1))]); its[its>itMax] = itMax
		lam = as.numeric(sapply(names(its), function(n){substring(strsplit(n, 'L')[[1]][2], 1, 4)}))
		W = as.numeric(sapply(names(its), function(n){substring(strsplit(n, 'W')[[1]][2], 1, 5)}))
		zMax = out[[i]]$Zmax[its]
		xMax = out[[i]]$Xmax[its,]
		colnames(xMax) = sprintf('x%d', 1:dm)
		#
		runsTil = c()
		for(ii in 1:length(its)){
			#
			if( is.na(its[ii]) ){ r=NA 
			}else{ r=its[ii]-min(which(out[[i]]$Zmax==zMax[ii])) }
			#
			runsTil = c(runsTil, r)
		}
		#
		nr = length(its)
		outter = rbind(outter, cbind(lam, W, its, zMax, xMax, runsTil))
	}
	#
	return( as.data.frame(outter) )
}

#
unpackAuto = function(out){
	#
	outter = c()#data.frame(lamda=numeric(), W=integer(), it=integer())
	ii = 1
	for(i in 1:length(out)){
		#
		nome = sprintf('ewmaAutoW%.2f', wGrid)
		#
		its = unlist(out[[i]]$itConv[nome]); its[its>itMax] = itMax
		#lam = as.numeric(sapply(as.character(wGrid), function(n){}))
		W = wGrid
		zMax = out[[i]]$Zmax[its]
		xMax = out[[i]]$Xmax[its,]
		colnames(xMax) = sprintf('x%d', 1:dm)
		#
                runsTil = c()
                for(ii in 1:length(its)){
                        #
                        if( is.na(its[ii]) ){ r=NA 
                        }else{ r=its[ii]-min(which(out[[i]]$Zmax==zMax[ii])) }
                        #
                        runsTil = c(runsTil, r)
                }
		#
		nr = length(its)
		outter = rbind(outter, cbind(W, its, zMax, xMax, runsTil))
	}
	#
	return( as.data.frame(outter) )
}

#
#
#

#
#load('testRun.RData')
#load('testL0.100.65W20.0070.00.RData')
#load('rastriginTestL0.100.65W20.0070.00.RData')
load('rastriginTestL0.100.65W20.00150.00.RData')
#load('rosenbrockTestL0.100.65W20.0070.00.RData')

#
rosenbrock = function(x){
        #
        x = matrix(x, ncol=2)
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

#
rastrigin = function(x, p=2){
        #
        x = matrix(x, ncol=p)
        out = matrix(NaN, nrow=dim(x)[1], ncol=1)
        for (i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = 10*length(ex)+sum(ex^2-10*cos(2*pi*ex))
        }
        return(out)
}
#
zThresh = 1e-4#1e-5
xThresh = rectVol*0.01^dm

#
#AUTO
#

#
for(w in wGrid){
	#
	len = length(out[[1]]$autoLams[[as.character(w)]])
	autoLamMat = matrix(NA, nrow=M, ncol=itMax-w)
	autoLamMat[1,1:len] = out[[1]]$autoLams[[as.character(w)]]
	#
	png(sprintf('%sAutoLam%s.png', name, w))
	plot(1:itMax, c(rep(NA, w), out[[1]]$autoLams[[as.character(w)]], rep(NA, (itMax-w)-len)), type='l', 
		col=adjustcolor("black", alpha.f = 0.5),
		ylim=c(0, 0.75),
		ylab="Lambda Estimate",
		main=sprintf("W Fixed at %s; Lambda Estimated", w)
	)
	#points(out[[1]]$itConv['ewmaAuto']-W, out[[1]]$autoLams[ unlist(out[[1]]$itConv['ewmaAuto']-W) ], col='red', pch=19)
	for(i in 2:M){ 
		#
		lam = out[[i]]$autoLams[[as.character(w)]]
		len = length(lam)
		autoLamMat[i,1:len] = lam
		#
		lines(1:itMax, c(rep(NA, w), lam, rep(NA, (itMax-w)-len)), col=adjustcolor("black", alpha.f = 0.5)) 
		#points(out[[i]]$itConv['ewmaAuto']-W, out[[i]]$autoLams[ unlist(out[[i]]$itConv['ewmaAuto']-W) ], col='red', pch=19)
	}
	#
	lines(1:itMax, c(rep(NA, w), colMeans(autoLamMat, na.rm=T)), col='blue', lwd=3)
	dev.off()
}

#
#
#

#
itX = unpackIts(out)
xNorm = apply(itX[,c('x1', 'x2')]-xMin, 1, Norm)
zDist = itX[,'zMax']-zMin
itX = cbind(itX, zDist, xNorm)

#
png(sprintf('%sItPairs.png', name))
pairs(itX[,-which(colnames(itX)%in%c('x1', 'x2', 'zMax'))])
dev.off()

##
#scatterplot3d(itX$W, itX$lam,  itX$its/(itMax+1))
#dev.new()
#scatterplot3d(itX$lam, itX$W,  itX$its/(itMax+1))

#
threshIt = sapply( out, function(x){ x$itConv[[1]] } )
threshZ = sapply( out, function(x){ x$Zmax[x$itConv[[1]]] } )
threshX = t(sapply( out, function(x){ x$Xmax[x$itConv[[1]],] } ))
#
threshBetaZ = mean(threshZ>zThresh)
threshBetaX = mean(threshX>xThresh)

#
#PROBS
#

#betaAvgZ = mean(itX[itX$lam==0.45 & itX$W==35,'zDist']>zThresh, na.rm=T)
FPR = mean(itX[itX$lam==0.45 & itX$W==35,'zDist']>zThresh, na.rm=T)
TP = sum(itX[itX$lam==0.45 & itX$W==35,'zDist']<zThresh)

#
betasZ = matrix(NA, length(wGrid), length(lamGrid))
colnames(betasZ) = sprintf('L:%.2f', lamGrid)
rownames(betasZ) = sprintf('W:%.2f', wGrid)
for(j in 1:ncol(betasZ)){
	for(i in 1:nrow(betasZ)){
		betasZ[i, j] = mean(itX[itX$lam==round(lamGrid[j],2) & itX$W==wGrid[i],'zDist']>zThresh, na.rm=T)
	}
}

#
betaAvgX = mean(itX[itX$lam==0.45 & itX$W==35,'xNorm']>xThresh, na.rm=T)
betasX = matrix(NA, length(wGrid), length(lamGrid))
colnames(betasX) = sprintf('L:%.2f', lamGrid)
rownames(betasX) = sprintf('W:%.2f', wGrid)
for(j in 1:ncol(betasX)){
	for(i in 1:nrow(betasX)){
		betasX[i, j] = mean(itX[itX$lam==round(lamGrid[j],2) & itX$W==wGrid[i],'xNorm']>xThresh, na.rm=T)
	}
}

#
png(sprintf('%sFaslePosZ.png', name))
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(t(betasZ), col=function(n) hcl.colors(n, "Blue-Red 3"))
dev.off()
#
png(sprintf('%sFaslePosX.png', name))
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(t(betasX), col=function(n) hcl.colors(n, "Blue-Red 3"))
dev.off()

##
#zProbW = sapply(wGrid, function(x){ mean(itX[itX$lam==0.4 & itX$W==x,'zDist']>zThresh, na.rm=T) })
#zProbL = sapply(lamGrid, function(x){ mean(itX[itX$lam==x & itX$W==40,'zDist']>zThresh, na.rm=T) })
##
#xProbW = sapply(wGrid, function(x){ mean(itX[itX$lam==0.4 & itX$W==x,'xNorm']>xThresh, na.rm=T) })
#xProbL = sapply(lamGrid, function(x){ mean(itX[itX$lam==x & itX$W==40,'xNorm']>xThresh, na.rm=T) })

#
itAuto = unpackAuto(out)
xNorm = apply(itAuto[,c('x1', 'x2')]-xMin, 1, Norm)
zDist = itAuto[,'zMax']-zMin
itAuto = cbind(itAuto, zDist, xNorm)
#
zProbAutoW = sapply(wGrid, function(x){ mean(itAuto[itAuto$W==x,'zDist']>zThresh, na.rm=T) })
xProbAutoW = sapply(wGrid, function(x){ mean(itAuto[itAuto$W==x,'xNorm']>xThresh, na.rm=T) })
#
png(sprintf('%sAutoFaslePosZThresh%.1e.png', name, zThresh))
plot(wGrid, zProbAutoW, 'l', lwd=3, ylim=c(0,1), main=sprintf('False Positive w/ Threshhold=%.2e', zThresh))
#
i = 1
cols = brewer.pal(length(lamGrid), "RdYlBu")
for(l in lamGrid){
	zProbW = sapply(wGrid, function(x){ mean(itX[itX$lam==l & itX$W==x,'zDist']>zThresh, na.rm=T) })
	lines(wGrid, zProbW, 'l', lwd=3, col=cols[i])
	i = i+1
}
lines(wGrid, zProbAutoW, 'l', lwd=5)
abline(h=threshBetaZ, lty=2, col='red')
dev.off()

#
#False negative
#

#Number of runs until convergence from best found before convergence (false negative)
alphasZ = matrix(NA, length(wGrid), length(lamGrid))
colnames(alphasZ) = sprintf('L:%.2f', lamGrid)
rownames(alphasZ) = sprintf('W:%.2f', wGrid)
for(j in 1:ncol(betasZ)){
        for(i in 1:nrow(betasZ)){
                alphasZ[i, j] = mean(itX[itX$lam==round(lamGrid[j],2) & itX$W==wGrid[i],'runsTil']/wGrid[i], na.rm=T)
        }
}

#
png(sprintf('%sFasleNegZ.png', name))
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(t(alphasZ), col=function(n) hcl.colors(n, "Blue-Red 3"))
dev.off()

#
png(sprintf('%sErrorZ.png', name))
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(t(alphasZ*betasZ), col=function(n) hcl.colors(n, "Blue-Red 3"))
dev.off()

#
#JUNK
#

##
#cOut = matrix(NA)
#for(i in 1:M){
#	#
#	zMax = out[[i]]$Zmax
#	xMax = out[[i]]$Xmax
#	#out[[i]]$Zmax[ itX[,'its'][i] ]
#	#out[[i]]$Xmax[ itX[,'its'][i], ]
#	#
#	for(w in wGrid){
#		#
#		for(l in lamGrid){
#			#
#			itty = unlist(out[[i]]$itConv[sprintf('ewmaL%.2fW%.2f', l, w)])
#			if(itty>itMax){ itty=itMax }
#			#itX[itX$lam==l & itX$W==w, 'its']
#			#if( zMax[itty]<zthresh ){
#				#	
#			#}
#			print(c(i, w, l))
#			print(itty)
#			XNorm = Norm( xMax[itty,]-xMin )
#			print(zMax[itty])
#			print(XNorm)
#			#if( XNorm<xThresh ){
#				#	
#			#}
#		}
#	}
#}







