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
unpackAuto = function(out){
        #
        outter = c()#data.frame(lamda=numeric(), W=integer(), it=integer())
        ii = 1
        for(i in 1:length(out)){
                #
                nome = sprintf('ewmaAutoW%.2f', wGrid)
                #
                its = unlist(out[[i]]$itConv[nome]);
                #
                #print( its )
                #print( as.matrix(out[[i]]$Xmax) )
                #its[its>itMax] = itMax
                #lam = as.numeric(sapply(as.character(wGrid), function(n){}))
                W = wGrid
                zCvg = out[[i]]$Zmax[its]
                xMax = matrix(matrix(out[[i]]$Xmax, ncol=dm)[its,], ncol=dm)
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
#LOAD
#

#
#load('rastriginNoCensorL0.1000.65W20.0080.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})>2; out=out[isGood]; M=sum(isGood); itMax=500
#load('rosenbrockNoCensorL0.1000.65W20.0040.00M100.RData')
#load('rosenbrockNoCensorL0.1000.65W10.0040.00M100.RData')
#load('rosenbrockNoCensorL0.1000.65W5.0045.00M100.RData')
load('./zooid3/rosenbrock3DL0.1000.65W5.0050.00M100.RData')
#load('./zooid1/grlee12L0.1000.65W20.0040.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})>2; out=out[isGood]; M=sum(isGood); itMax=500
#load('./zooid1/grlee12L0.1000.65W10.0040.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})>2; out=out[isGood]; M=sum(isGood); itMax=500
#load('./zooid1/grlee12L0.1000.65W5.0040.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})>2; out=out[isGood]; M=sum(isGood); itMax=500
#load('./zooid2/michal2DL0.1000.65W20.0040.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})!=0; out=out[isGood]; M=sum(isGood); itMax=500
#load('./zooid2/michal2DL0.1000.65W10.0040.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})!=0; out=out[isGood]; M=sum(isGood); itMax=500
#load('./zooid2/michal5DL0.1000.65W20.0040.00M100.RData')
#load('./zooid2/michal5DL0.1000.65W50.00100.00M100.RData')
#load('./zooid2/michal3DL0.1000.65W100.00200.00M100.RData'); isGood=sapply(out, function(x){length(names(x))})>2; out=out[isGood]; M=sum(isGood); itMax=500
#load('./zooid4/levyL0.1000.65W20.0040.00.RData')
#load('./zooid4/levyL0.1000.65W40.0080.00.RData')
#load('./zooid4/levy2DL0.1000.65W5.0099.00.RData')

#
#FRONT MATTER
#

#truth threshhold
zThresh = 1e-4

#
#EI
#

#
threshs = sprintf('thresh%1.1e', thresholdGrid)
threshZ = t(sapply( out, function(x){ x$Zmax[unlist(x$itConv[threshs])] } ))
colnames(threshZ) = threshs
threshFPR = colMeans((threshZ-zMin[1])>zThresh)
threshFPRSE = sqrt(apply((threshZ-zMin[1])>zThresh, 2, var)/nrow(threshZ))

#
ts = c(thresholdGrid, rev(thresholdGrid))

#FPR
png(sprintf('%sEIFPRZThresh%.1e.png', name, zThresh))
plot(thresholdGrid, threshFPR, pch=20)
polygon(ts, c(threshFPR+2*threshFPRSE, rev(threshFPR-2*threshFPRSE)),
        col=adjustcolor("black", alpha.f=0.2),
        border=F
)
dev.off()
#lines(thresholdGrid, arlCThreshTil, 'l', lwd=3)

#
threshIt = t(sapply( out, function(x){ x$itConv[threshs] } ))
#
threshRunsTils = c()
#
arlCThreshTil = c()
arlCThreshTilSE = c()
#
arlCThreshTilGiven = c()
arlCThreshTilGivenSE = c()
#
arlCThreshTilChop = c()
arlCThreshTilChopSE = c()
for(thresh in threshs){
        threshRunsTil = c()
        for(ii in 1:length(threshIt[,thresh])){
                #
                less = out[[ii]]$Zmax<=zThresh
                if( is.na(threshIt[ii,thresh]) | !any(less) ){ r=NA
                }else{ r=unlist(threshIt[ii,thresh])-min(which(less)) }
                #NOTE: negative arl means optimization entered the truth region was found after alg claimed convergence
                threshRunsTil = c(threshRunsTil, r)
        }
        #
        threshRunsTils = cbind(threshRunsTils, threshRunsTil)
        #
	cond = threshZ<zThresh
	arlCThreshTil = c(arlCThreshTil, mean(threshRunsTil[cond], na.rm=T))
        arlCThreshTilNum = sum(!is.na(threshRunsTil[cond]))
        arlCThreshTilSE = c(arlCThreshTilSE, sqrt( var(threshRunsTil[cond], na.rm=T)/arlCThreshTilNum ))
	#E[ runsTil*1(runsTil>0) ]: use same cond from above
	threshRunsTilChop = threshRunsTil
	threshRunsTilChop[threshRunsTilChop<0] = 0 
	arlCThreshTilChop = c(arlCThreshTilChop, mean(threshRunsTilChop[cond], na.rm=T))
        arlCThreshTilChopNum = sum(!is.na(threshRunsTilChop[cond]))
        arlCThreshTilChopSE = c(arlCThreshTilChopSE, sqrt( var(threshRunsTilChop[cond], na.rm=T)/arlCThreshTilChopNum ))
	#E[runsTil | runsTill>0]
	cond = threshZ<zThresh & threshRunsTil>=0
	arlCThreshTilGiven = c(arlCThreshTilGiven, mean(threshRunsTil[cond], na.rm=T))
        arlCThreshTilGivenNum = sum(!is.na(threshRunsTil[cond]))
        arlCThreshTilGivenSE = c(arlCThreshTilGivenSE, sqrt( var(threshRunsTil[cond], na.rm=T)/arlCThreshTilGivenNum ))	
}
colnames(threshRunsTils) = threshs
#
names(arlCThreshTil) = threshs
names(arlCThreshTilSE) = threshs
#
names(arlCThreshTilChop) = threshs
names(arlCThreshTilChopSE) = threshs
#
names(arlCThreshTilGiven) = threshs
names(arlCThreshTilGivenSE) = threshs

#
png(sprintf('%sEIARLZThresh%.1e.png', name, zThresh))
plot(thresholdGrid, arlCThreshTil, pch=20, ylim=c(-40, 120))#, 'l')
polygon(ts, c(arlCThreshTil+2*arlCThreshTilSE, rev(arlCThreshTil-2*arlCThreshTilSE)),
        col=adjustcolor("black", alpha.f=0.2),
        border=F
)
points(thresholdGrid, arlCThreshTilChop, pch=20, col='red')
polygon(ts, c(arlCThreshTilChop+2*arlCThreshTilChopSE, rev(arlCThreshTilChop-2*arlCThreshTilChopSE)),
        col=adjustcolor("red", alpha.f=0.2),
        border=F
)
points(thresholdGrid, arlCThreshTilGiven, pch=20, col='blue')
bounds = c(arlCThreshTilGiven+2*arlCThreshTilGivenSE, rev(arlCThreshTilGiven-2*arlCThreshTilGivenSE))
notBroke = !is.na(bounds)
polygon(ts[notBroke], bounds[notBroke],
        col=adjustcolor("blue", alpha.f=0.2),
        border=F
)
legend('topright', legend=c('E[RL]', 'E[RL x 1(RL>0)]', 'E[RL | RL>=0]'), col=c('black', 'red', 'blue'), pch=20)
dev.off()
#lines(thresholdGrid, arlCThreshTil, 'l', lwd=3)

#
#EWMA ELAI
#

#
itAuto = unpackAuto(out)
zDist = itAuto[,'zCvg']-zMin
itAuto = cbind(itAuto, zDist)

# & itAuto$its<itMax
#we want arlNC as large as possible
#in the case of no convergence, and the run gets to itMax, the run did the correct thing ans the true its-value may actually be larger
#in the case of no convergence removing itAuto$its>=itMax makes the ARL articifically low which gives a generously conservative of the ARL 
arlNC = sapply(wGrid, function(x){ mean(itAuto$its[itAuto$W==x & (itAuto$zCvg-zMin[1])>zThresh], na.rm=T) })
arlNCNum = sapply(wGrid, function(x){ sum(itAuto$its[itAuto$W==x & (itAuto$zCvg-zMin[1])>zThresh]>=0, na.rm=T) })
arlNCSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$its[itAuto$W==x & (itAuto$zCvg-zMin[1])<zThresh], na.rm=T)) })
arlNCSE  = arlNCSE/sqrt(arlNCNum)
# & itAuto$its<itMax
#we want arlCTil as small as possible
#in the case of convergence, and the run gets to itMax, the run length is likely not doing the correct thing, and may have been larger than itMax
#if itAuto$its<itMax limitation is not included the right-censored itMax its-value biases the arlCTil to be smaller than it should be (exactly the wrong direction)
cond = (itAuto$zCvg-zMin[1])<zThresh
arlCTil  = sapply(wGrid, function(x){ mean(itAuto$runsTil[itAuto$W==x & cond], na.rm=T) })
arlCNumTil = sapply(wGrid, function(x){ sum(itAuto$runsTil[itAuto$W==x & cond]>=0, na.rm=T) })
arlCTilSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$runsTil[itAuto$W==x & cond], na.rm=T)) })
arlCTilSE  = arlCTilSE/sqrt(arlCNumTil)
#E[ runsTil*1(runsTil>0) ]
itAuto$runsTilChop = itAuto$runsTil
itAuto$runsTilChop[itAuto$runsTilChop<0] = 0
arlCTilChop  = sapply(wGrid, function(x){ mean(itAuto$runsTilChop[itAuto$W==x & cond], na.rm=T) })
arlCNumTilChop = sapply(wGrid, function(x){ sum(itAuto$runsTilChop[itAuto$W==x & cond]>=0, na.rm=T) })
arlCTilChopSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$runsTilChop[itAuto$W==x & cond], na.rm=T)) })
arlCTilChopSE  = arlCTilChopSE/sqrt(arlCNumTilChop)
#E[runsTil | runsTill>0]
cond = (itAuto$zCvg-zMin[1])<zThresh & itAuto$runsTil>=0
arlCTilGiven  = sapply(wGrid, function(x){ mean(itAuto$runsTil[itAuto$W==x & cond], na.rm=T) })
arlCNumTilGiven = sapply(wGrid, function(x){ sum(itAuto$runsTil[itAuto$W==x & cond]>=0, na.rm=T) })
arlCTilGivenSE  = sapply(wGrid, function(x){ sqrt(var(itAuto$runsTil[itAuto$W==x & cond], na.rm=T)) })
arlCTilGivenSE  = arlCTilGivenSE/sqrt(arlCNumTilGiven)

#
dubs = c(wGrid, rev(wGrid))
int = c(arlCTil+2*arlCTilSE, rev(arlCTil-2*arlCTilSE))
show = !is.na(int)
#int = c(arlNC+2*arlNCSE, rev(arlNC-2*arlNCSE))
#show = !is.na(int)

#ARL
png(sprintf('%sELAIARLZThresh%.1e.png', name, zThresh))
plot(wGrid, arlCTil,
        lwd=3,
        type='l',
        ylim=range(c(arlCThreshTil, int), na.rm=T), #c(min(arlCThreshTil), max(int[show])), 
        #axes=F,
        xlab='',
        ylab='',
)
polygon(dubs[show], int[show],
        col=adjustcolor("black", alpha.f=0.2),
        border=F
)
lines(wGrid, arlCTil, lwd=3)
#
lines(wGrid, arlCTilChop, lwd=3, col='red')
intChop = c(arlCTilChop+2*arlCTilChopSE, rev(arlCTilChop-2*arlCTilChopSE))
showChop = !is.na(intChop)
polygon(dubs[showChop], intChop[showChop],
        col=adjustcolor("red", alpha.f=0.2),
        border=F
)
#
lines(wGrid, arlCTilGiven, lwd=3, col='blue')
intGiven = c(arlCTilGiven+2*arlCTilGivenSE, rev(arlCTilGiven-2*arlCTilGivenSE))
showGiven = !is.na(intGiven)
polygon(dubs[showGiven], intChop[showGiven],
        col=adjustcolor("blue", alpha.f=0.2),
        border=F
)
legend('bottomright', legend=c('E[RL]', 'E[RL x 1(RL>0)]', 'E[RL | RL>=0]'), col=c('black', 'red', 'blue'), lwd=3)
#axis(side=4, at=pretty(range(c(arlCThreshTil, int))))
#mtext('ARL', side=4, line=3)
dev.off()

#formerly zProbAutoW, and zProbAutoWSD
ewmaFPR = sapply(wGrid, function(x){ mean(itAuto[itAuto$W==x,'zDist']>zThresh, na.rm=T) })
ewmaFPRSE = sapply(wGrid, function(x){ sqrt(var(itAuto[itAuto$W==x,'zDist']>zThresh, na.rm=T)/sum(!is.na(itAuto[itAuto$W==x,'zDist']>zThresh))) })

#FPR
png(sprintf('%sELAIFPRZThresh%.1e.png', name, zThresh))
plot(wGrid, ewmaFPR, 'l',
        lwd=3,
        ylim=c(0,1),
        main=sprintf('Errors w/ Truth Threshhold=%.2e', zThresh),
        ylab='FPR'
)
polygon(dubs, c(ewmaFPR+2*ewmaFPRSE, rev(ewmaFPR-2*ewmaFPRSE)),
        col=adjustcolor("black", alpha.f=0.2),
        border=F
)
lines(wGrid, ewmaFPR, 'l', lwd=3)
dev.off()


#
#rowCombine FRP
#

#
cols = brewer.pal(8, "Set1")

#
png(sprintf('%sFPRCombine%.1e.png', name, zThresh))
par(mar=c(5, 4, 5, 4)+0.3)

#EWMA
plot(wGrid, ewmaFPR, 'l',
        lwd=3,
        ylim=c(0,1),
        main='', #sprintf('Errors w/ Truth Threshhold=%.2e', zThresh),
	#axes=F,
	xlab='',
        ylab='FPR',
	col=cols[1]
)
polygon(dubs, c(ewmaFPR+2*ewmaFPRSE, rev(ewmaFPR-2*ewmaFPRSE)),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
lines(wGrid, ewmaFPR, 'l', lwd=3, col=cols[1])

#
mtext('Window Size', side=1, line=3, col=cols[1])

#
par(new=T)

#THRESH
plot(-thresholdGrid, threshFPR, 'l', 
	ylab="", 
	xlab="", 
	axes=F, 
	col=cols[2], 
	lwd=3,
	ylim=c(0,1)
)
polygon(-ts, c(threshFPR+2*threshFPRSE, rev(threshFPR-2*threshFPRSE)),
        col=adjustcolor(cols[2], alpha.f=0.2),
        border=F
)

#
ax2 = pretty(range(c(-thresholdGrid, -ts)))
axis(side=3, at=ax2, labels=-ax2)
mtext('EI Threshold', side=3, line=3, col=cols[2])

#
dev.off()

#
#rowCombine ARL
#

#
png(sprintf('%sARLCombine%.1e.png', name, zThresh))
par(mar=c(5, 4, 5, 4)+0.3)

#ewma interval
intChop = c(arlCTilChop+2*arlCTilChopSE, rev(arlCTilChop-2*arlCTilChopSE))
showChop = !is.na(intChop)
#thresh interval
thInt = c(arlCThreshTilChop+2*arlCThreshTilChopSE, rev(arlCThreshTilChop-2*arlCThreshTilChopSE))
ylim = range(intChop, thInt, na.rm=T)
#
plot(wGrid, arlCTilChop, 'l',
	lwd=3, 
	col=cols[1],
	xlab='',
	ylab="ARL",
	ylim=ylim
)
polygon(dubs[showChop], intChop[showChop],
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)

#
mtext('Window Size', side=1, line=3, col=cols[1])


#
par(new=T)

#THRESH
plot(-thresholdGrid, arlCThreshTilChop, 'l', 
	col=cols[2],
	lwd=3,
	xlab='',
	ylab='',
	axes=F,
	ylim=ylim
)
polygon(-ts, thInt,
        col=adjustcolor(cols[2], alpha.f=0.2),
        border=F
)

#
ax2 = pretty(range(c(-thresholdGrid, -ts)))
axis(side=3, at=ax2, labels=-ax2)
mtext('EI Threshold', side=3, line=3, col=cols[2])

#
dev.off()



