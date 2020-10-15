rm(list=ls())

#
#FUNCTIONS
#

getN = function(lm, lv){
	#lm: the lognormal mean
	#lv: the lognormal variance
	return( exp(lm+(lv/2)) )
}

#
#MAIN
#

##
##load('rosenbrockThirtyMPISim.RData')
##load('rosenbrockThirtyMPISimW10.RData')
##load('rosenbrockThirtyMPISimW20.RData')
#load('rosenbrockThirtyMPISimW30.RData')
#
#load('rastriginWFiftyTestNNSim.RData')
load('rastriginW90NNSimW90.RData')
#
itConv = c()
ZBothMins = c()
ZEwmaMins = c()
ZThreshMins = c()
#
skp = c()
#Xorigs = c()
#Zorigs = c()
#Xtrans = c()
#Ztrans = c()
#Xmins = c()
##
origThresh = 5e-4
transThresh = -9.5 #6e-4#1e-4
##
#dev.new()
#plot(clusterOut[[1]]$origMean, ylim=c(0, 0.5), xlim=c(0, 350))
#dev.new()
#plot(clusterOut[[1]]$transMean, xlim=c(0, 350), ylim=c(-16, 0) )
for(i in 1:M){
	#broken optimizer
	if( typeof(clusterOut[[i]])=="character" ){ 
		#
		skp = c(skp, i)
		next
	}
	#
	#origWhere = which(clusterOut[[i]]$origMean<origThresh)[1]+40-1
	#transWhere = which(clusterOut[[i]]$transMean<transThresh)[1]+40-1
	##print(sprintf("%f %f", origWhere, transWhere))
	##Zorigs = 
	##Ztrans =
	#if(!is.na(origWhere)){ 
	#	Xorigs=rbind(Xorigs, clusterOut[[i]]$X[origWhere,])
	#	Zorigs=c(Zorigs, clusterOut[[i]]$Zmax)
	#}
	#if(!is.na(transWhere)){ 
	#	Xtrans=rbind(Xtrans, clusterOut[[i]]$X[transWhere,]) 
	#	#Ztrans=rbind()
	#}
	#
	its = clusterOut[[i]]$itConv
	itConv = rbind(itConv, its)
	#
	ZBothMin = min(clusterOut[[i]]$Zmax)
	ZBothMins = c(ZBothMins, ZBothMin)
	#
	ZEwmaMin = min(clusterOut[[i]]$Zmax[1:its[1]])
	ZEwmaMins = c(ZEwmaMins, ZEwmaMin)
	#
	ZThreshMin = min(clusterOut[[i]]$Zmax[1:its[2]])
	ZThreshMins = c(ZThreshMins, ZThreshMin)
	#
	#Xwhere = which(clusterOut[[i]]$Zmax==Zmin)[1]+40-2 #40 is the size of the initial set, one index over for the initial set and one for actual optimization
	#Xmins = rbind(Xmins, clusterOut[[i]]$X[Xwhere,])	
	##	
	##
	#points(1:length(clusterOut[[i]]$origMean), clusterOut[[i]]$origMean)
	#points(1:length(clusterOut[[i]]$transMean), clusterOut[[i]]$transMean)
}
#
#pdf('itConv.pdf')
pdf('itConvRast.pdf')
plot(itConv,
	xlim=c(0,300), 
	ylim=c(0,300),
	main='Rosenbrock' # 'Rastrigin'
)
lines(-100:400, -100:400)
dev.off()
#
bothCol = max.col(itConv)
d = dim(itConv)[1]
bothConv = matrix(NaN, nrow=d, ncol=1)
for(i in 1:d){ bothConv[i]=itConv[i,bothCol[i]] }
#dev.new()
#plot(itConv[,1], bothConv)
#dev.new()
#plot(itConv[,2], bothConv)
#
Zthresh = 1e-4 #0.1#
Neff = length(ZBothMins)
#
ZBothWrong = sum(ZBothMins>Zthresh)/Neff
ZEwmaWrong = sum(ZEwmaMins>Zthresh)/Neff
ZThreshWrong = sum(ZThreshMins>Zthresh)/Neff
#
#pdf('ewmaSimHist.pdf')
pdf('ewmaSimHistRast.pdf')
hist(ZEwmaMins, 100,
	main=sprintf('\nMean: %0.6g\nMedian: %0.6g\nVar: %0.6g\nPr(Min Z>%0.6g) = %0.6g', mean(ZEwmaMins), median(ZEwmaMins), var(ZEwmaMins), Zthresh, ZEwmaWrong),
	xlab='EWMA Minimum Z'#,
	#xlim=c(0,0.0025)
	#xlim=c(0,)
)
dev.off()
#
#pdf('threshSimHist.pdf')
pdf('threshSimHistRast.pdf')
hist(ZThreshMins, 100,
	main=sprintf('\nMean: %0.6g\nMedian: %0.6g\nVar: %0.6g\nPr(Min Z>%0.6g) = %0.6g', mean(ZThreshMins), median(ZThreshMins), var(ZThreshMins), Zthresh, ZThreshWrong),
	xlab='Threshold Minimum Z'#,
	#xlim=c(0,0.0025)
	#xlim=c(0,)
)
dev.off()
#
#pdf('bothSimHist.pdf')
pdf('bothSimHistRast.pdf')
hist(ZBothMins, 100,
	main=sprintf('\nMean: %0.6g\nMedian: %0.6g\nVar: %0.6g\nPr(Min Z>%0.6g) = %0.6g', mean(ZBothMins), median(ZBothMins), var(ZBothMins), Zthresh, ZBothWrong),
	xlab='Both Minimum Z',
	#xlim=c(0,0.0025)
	#xlim=c(0,)
)
dev.off()




#print(sprintf('Both: %f', ZBothWrong))
#print(sprintf('EWMA: %f', ZEwmaWrong))
#print(sprintf('Thresh %f', ZThreshWrong))
##
#Xthresh = 1
#Xnorms = apply(Xmins, 1, function(x){sqrt(sum(x^2))})
#Xwrong = sum(Xnorms>Xthresh)/Neff
##
#Xon = apply(Xorigs, 1, function(x){sqrt(sum(x^2))})
##Xtn = apply(Xtrans, 1, function(x){sqrt(sum(x^2))})
#
##Zmins = unlist(lapply(clusterOut, function(x){ min(x[[3]]) }))


