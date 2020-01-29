rm(list=ls())

#
#
#

#
load('./bh_1800_b0_3000_M2_040_300RDat.RData')
par = read.csv('bh_1800_b03000_M2_040_300Par.csv')
#
sSum = matrix(NA, nrow=length(outs), ncol=4)
colnames(sSum) = c('lo', 'hi', 'median', 'mean')
qSum = matrix(NA, nrow=length(outs), ncol=4)
colnames(qSum) = c('lo', 'hi', 'median', 'mean')
MSum = matrix(NA, nrow=length(outs), ncol=4)
colnames(MSum) = c('lo', 'hi', 'median', 'mean')
deltaSum = matrix(NA, nrow=length(outs), ncol=4)
colnames(deltaSum) = c('lo', 'hi', 'median', 'mean')
i = 1
for(o in outs){
	#
	sam = extract(o, permuted=T)
	#
	deltaSum[i,c('lo', 'median', 'hi')] = quantile(sam$delta, probs=c(0.025, 0.5, 0.975))
	deltaSum[i,'mean'] = mean(sam$delta)
	#
	MSum[i,c('lo', 'median', 'hi')] = quantile(-log(1-sam$delta), probs=c(0.025, 0.5, 0.975))
	MSum[i,'mean'] = mean(-log(1-sam$delta))
	#
	qSum[i,c('lo', 'median', 'hi')] = quantile(sam$q, probs=c(0.025, 0.5, 0.975))
	qSum[i,'mean'] = mean(sam$q)
	#
	sSum[i,c('lo', 'median', 'hi')] = quantile(sam$s, probs=c(0.025, 0.5, 0.975))
	sSum[i,'mean'] = mean(sam$s)
	#
	i = i + 1
}
#
M = seq(0.04, 0.3, 0.005)
##
#plot(M, deltaSum[,'mean'], 'l', 
#	lwd=3,
#	ylim=c(min(deltaSum), max(deltaSum))
#)
#polygon(c(M, rev(M)), c(deltaSum[,'hi'], rev(deltaSum[,'lo'])))
##
dev.new()
plot(M, MSum[,'mean'], 'l', 
	lwd=3,
	ylim=c(min(MSum), max(MSum))
)
polygon(c(M, rev(M)), c(MSum[,'hi'], rev(MSum[,'lo'])))
#
dev.new()
plot(M, qSum[,'mean'], 'l', 
	lwd=3,
	ylim=c(min(qSum), max(qSum))
)
polygon(c(M, rev(M)), c(qSum[,'hi'], rev(qSum[,'lo'])))
#
dev.new()
plot(M, sSum[,'mean'], 'l', 
	lwd=3,
	ylim=c(min(sSum), max(sSum))
)
polygon(c(M, rev(M)), c(sSum[,'hi'], rev(sSum[,'lo'])))






##
#glob0 = Sys.glob('pics/bh_1800_b0_3000_M2_040_300RDat\ *.RData')
#glob1 = Sys.glob('pics/bh_1800_b0_3000_M2_040_300RDat1*.RData')
#glob2 = Sys.glob('pics/bh_1800_b0_3000_M2_040_300RDat2*.RData')
#glob3 = Sys.glob('pics/bh_1800_b0_3000_M2_040_300RDat3*.RData')
##
#glob = c(glob0, glob1, glob2, glob3)
##
#
#for(g in glob){
#	#
#	load(g)
#}

