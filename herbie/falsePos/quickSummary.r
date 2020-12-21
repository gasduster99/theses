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
		zCvg = out[[i]]$Zmax[its]
		xMax = out[[i]]$Xmax[its,]
		colnames(xMax) = sprintf('x%d', 1:dm)
		#
		runsTil = c()
		for(ii in 1:length(its)){
			#
			less = out[[i]]$Zmax<=zThresh
			if( is.na(its[ii]) | !any(less)){ r=NA 
			}else{ r=its[ii]-min(which(less)) }
			#
			runsTil = c(runsTil, r)
		}
		#
		nr = length(its)
		outter = rbind(outter, cbind(lam, W, its, zCvg, xMax, runsTil))
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
		zCvg = out[[i]]$Zmax[its]
		xMax = out[[i]]$Xmax[its,]
		colnames(xMax) = sprintf('x%d', 1:dm)
		#
                runsTil = c()
                for(ii in 1:length(its)){
                        #
			less = out[[i]]$Zmax<=zThresh
                        if( is.na(its[ii]) | !any(less) ){ r=NA 
                        }else{ r=its[ii]-min(which(less)) }
                        #
                        runsTil = c(runsTil, r)
                }
		#
		nr = length(its)
		outter = rbind(outter, cbind(W, its, zCvg, xMax, runsTil))
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

##
##AUTO
##
#
##
#for(w in wGrid){
#	#
#	len = length(out[[1]]$autoLams[[as.character(w)]])
#	autoLamMat = matrix(NA, nrow=M, ncol=itMax-w)
#	autoLamMat[1,1:len] = out[[1]]$autoLams[[as.character(w)]]
#	#
#	png(sprintf('%sAutoLam%s.png', name, w))
#	plot(1:itMax, c(rep(NA, w), out[[1]]$autoLams[[as.character(w)]], rep(NA, (itMax-w)-len)), type='l', 
#		col=adjustcolor("black", alpha.f = 0.5),
#		ylim=c(0, 0.75),
#		ylab="Lambda Estimate",
#		main=sprintf("W Fixed at %s; Lambda Estimated", w)
#	)
#	#points(out[[1]]$itConv['ewmaAuto']-W, out[[1]]$autoLams[ unlist(out[[1]]$itConv['ewmaAuto']-W) ], col='red', pch=19)
#	for(i in 2:M){ 
#		#
#		lam = out[[i]]$autoLams[[as.character(w)]]
#		len = length(lam)
#		autoLamMat[i,1:len] = lam
#		#
#		lines(1:itMax, c(rep(NA, w), lam, rep(NA, (itMax-w)-len)), col=adjustcolor("black", alpha.f = 0.5)) 
#		#points(out[[i]]$itConv['ewmaAuto']-W, out[[i]]$autoLams[ unlist(out[[i]]$itConv['ewmaAuto']-W) ], col='red', pch=19)
#	}
#	#
#	lines(1:itMax, c(rep(NA, w), colMeans(autoLamMat, na.rm=T)), col='blue', lwd=3)
#	dev.off()
#}

#
#
#

#
itX = unpackIts(out)
xNorm = apply(itX[,c('x1', 'x2')]-xMin, 1, Norm)
zDist = itX[,'zCvg']-zMin
itX = cbind(itX, zDist, xNorm)
#itMax = max(itX$its)+100

#
png(sprintf('%sItPairs.png', name))
pairs(itX[,-which(colnames(itX)%in%c('x1', 'x2', 'zCvg'))])
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
threshRunsTil = c()
for(ii in 1:length(threshIt)){
        #
	less = out[[ii]]$Zmax<=zThresh
        if( is.na(threshIt[ii]) | !any(less) ){ r=NA
        }else{ r=threshIt[ii]-min(which(less)) }
	#
        threshRunsTil = c(threshRunsTil, r)
}
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
zDist = itAuto[,'zCvg']-zMin
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
#ARL
#

#
arlC  = sapply(wGrid, function(x){ mean(itAuto$its[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax], na.rm=T) })
arlCNum = sapply(wGrid, function(x){ sum(itAuto$its[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax]>=0, na.rm=T) })
arlCSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$its[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax], na.rm=T)) })
arlCSE  = arlCSE/sqrt(arlCNum)
#
arlNCTil = sapply(wGrid, function(x){ mean(itAuto$runsTil[itAuto$W==x & itAuto$zCvg>zThresh & itAuto$its<itMax], na.rm=T) })
arlNCNumTil = sapply(wGrid, function(x){ sum(itAuto$runsTil[itAuto$W==x & itAuto$zCvg>zThresh & itAuto$its<itMax]>=0, na.rm=T) })
arlNCTilSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$runsTil[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax], na.rm=T)) })
arlNCTilSE  = arlNCTilSE/sqrt(arlNCNumTil)

#the situation where its==itMax represents right censored observations

#we want arlNC as large as possible
#in the case of no convergence, and the run gets to itMax, the run did the correct thing ans the true its-value may actually be larger
#in the case of no convergence removing itAuto$its>=itMax makes the ARL articifically low which gives a generously conservative of the ARL 
arlNC = sapply(wGrid, function(x){ mean(itAuto$its[itAuto$W==x & itAuto$zCvg>zThresh & itAuto$its<itMax], na.rm=T) })
arlNCNum = sapply(wGrid, function(x){ sum(itAuto$its[itAuto$W==x & itAuto$zCvg>zThresh & itAuto$its<itMax]>=0, na.rm=T) })
arlNCSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$its[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax], na.rm=T)) })
arlNCSE  = arlNCSE/sqrt(arlNCNum)
#we want arlCTil as small as possible
#in the case of convergence, and the run gets to itMax, the run length is likely not doing the correct thing, and may have been larger than itMax
#if itAuto$its<itMax limitation is not included the right-censored itMax its-value biases the arlCTil to be smaller than it should be (exactly the wrong direction)
arlCTil  = sapply(wGrid, function(x){ mean(itAuto$runsTil[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax], na.rm=T) })
arlCNumTil = sapply(wGrid, function(x){ sum(itAuto$runsTil[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax]>=0, na.rm=T) })
arlCTilSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$runsTil[itAuto$W==x & itAuto$zCvg<zThresh & itAuto$its<itMax], na.rm=T)) })
arlCTilSE  = arlCTilSE/sqrt(arlCNumTil)

#
arlNCThresh = mean(threshIt[threshZ>zThresh])
arlNCThreshNum = length(threshIt[threshZ>zThresh])
arlNCThreshSE = sqrt( var(threshIt[threshZ>zThresh])/arlNCThreshNum )
#
arlCThreshTil = mean(threshRunsTil[threshZ<zThresh])
arlCThreshTilNum = length(threshRunsTil[threshZ<zThresh])
arlCThreshTilSE = sqrt( var(threshRunsTil[threshZ>zThresh])/arlCThreshTilNum )

#mean(itAuto$its[itAuto$W==x & itAuto$zCvg>zThresh & itAuto$its<itMax], na.rm=T)

#
falsePos = 1/arlNC
falseNeg = 1/(1-arlC)
#
falsePosTil = 1/arlNCTil
falseNegTil = 1/(1-arlCTil)

##
#png(sprintf('arl%s.png', name))
#plot(wGrid, arlC, lwd=3, type='l', ylim=c(30, itMax), main='Average Run Length')
#dubs = c(wGrid, rev(wGrid))
#int = c(arlNC+2*arlNCSE, rev(arlNC-2*arlNCSE))
#show = !is.na(int)
#polygon(dubs[show], int[show], col='pink', border=F)
#polygon(dubs, c(arlC+2*arlCSE, rev(arlC-2*arlCSE)), col='grey', border=F)
#lines(wGrid, arlC, lwd=3)
#lines(wGrid, arlNC, lwd=3, col='red')
#legend('topright', legend=c('True Claim', 'False Claim'), lty=1, lwd=3, col=c('black', 'red'))
#dev.off()

#
png(sprintf('arlNum%s.png', name))
plot(wGrid, arlCNum, lwd=3, type='l', ylim=c(0, 100), main='Average Sample Size')
lines(wGrid, arlNCNum, lwd=3, col='red')
legend('topright', legend=c('True Claim', 'False Claim'), lty=1, lwd=3, col=c('black', 'red'))
dev.off()

##
#png(sprintf('arlTil%s.png', name))
#plot(wGrid, arlCTil, lwd=3, type='l', ylim=c(0, 100), main='Average Run Length')
#dubs = c(wGrid, rev(wGrid))
#int = c(arlNCTil+2*arlNCTilSE, rev(arlNCTil-2*arlNCTilSE))
#show = !is.na(int)
#polygon(dubs[show], int[show], col='pink', border=F)
#polygon(dubs, c(arlCTil+2*arlCTilSE, rev(arlCTil-2*arlCTilSE)), col='grey', border=F)
#lines(wGrid, arlCTil, lwd=3)
#lines(wGrid, arlNCTil, lwd=3, col='red')
#legend('topright', legend=c('True Claim', 'False Claim'), lty=1, lwd=3, col=c('black', 'red'))
#dev.off()

#
png(sprintf('arlNumTil%s.png', name))
plot(wGrid, arlCNumTil, lwd=3, type='l', ylim=c(0, 100), main='Average Sample Size')
lines(wGrid, arlNCNumTil, lwd=3, col='red')
legend('topright', legend=c('True Claim', 'False Claim'), lty=1, lwd=3, col=c('black', 'red'))
dev.off()

#
png(sprintf('arlCombine%s.png', name))
plot(wGrid, arlCTil, lwd=3, type='l', ylim=c(0, itMax), main='Average Run Length')
dubs = c(wGrid, rev(wGrid))
int = c(arlNC+2*arlNCSE, rev(arlNC-2*arlNCSE))
show = !is.na(int)
polygon(dubs[show], int[show], col='pink', border=F)
int = c(arlCTil+2*arlCTilSE, rev(arlCTil-2*arlCTilSE))
show = !is.na(int)
polygon(dubs[show], int[show], col='grey', border=F)
lines(wGrid, arlCTil, lwd=3)
lines(wGrid, arlNC, lwd=3, col='red')
abline(h=arlNCThresh, col='red')
abline(h=arlNCThresh-2*arlNCThreshSE, col='red', lty=2)
abline(h=arlNCThresh+2*arlNCThreshSE, col='red', lty=2)
abline(h=arlCThreshTil)
abline(h=arlCThreshTil-2*arlCThreshTilSE, lty=2)
abline(h=arlCThreshTil+2*arlCThreshTilSE, lty=2)
legend('topright', legend=c('True Claim', 'False Claim'), lty=1, lwd=3, col=c('black', 'red'))
dev.off()

#itX$its[itX$zCvg<zThresh & itX$its<itMax]

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







