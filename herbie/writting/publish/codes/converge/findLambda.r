rm(list=ls())

source('spcSeeEWMANormFunc.r')

##name = "RosegpPic"; uu=30; tMin=0; 	#lambda=0.2
#name = "RoseEasyEasy"; uu=30; tMin=0;	#it=74		#lambda=0.2
name = "RastHard"; uu=40; tMin=0;	#it=96		#lambda=0.3
#
load(sprintf("seaSave%s.RData", name))
W = uu
it = 96
h = 0.01
#
ewmaOut = findEWMA(it, h)
ssOut = findSS(it, h)
#
ewmaPlot(ewmaOut, it)
dev.new()
plot(ssOut[[1]], ssOut[[2]])
abline(v=ewmaOut[[10]], col='red')








#ls = seq(0.01, 1, 0.01)
#nl = length(ls)
#ils = seq(1, nl)
#ss = matrix(NA, nrow=nl, ncol=1)
##
#out = ewmaConvChart(lambda, 50)
#ssError(out, lambda)
##
#for(i in ils){
#	ewma = ewmaConvChart(ls[i], it)
#	ss[i] = ssError(ewma, ls[i])
#	if(ss[i]<=min(ss[1:i])){ ewmaOut=ewma; lOut=ls[i] }
#}
