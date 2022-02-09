rm(list=ls())

#
library(deming)

#
#FUNCTIONS
#

#
pCube = function(t1, t2){ (t1-t2)^3*(t1>t2) }
pDiff = function(t1, t2){ (t1-t2)*(t1>t2) }
pInd = function(t1, t2){as.numeric(t1>t2)}
pInt = function(ti, tj){ ti * (ti/2-tj) * pInd(ti, tj) }


#
#MAIN
#

#
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(catch)
time = 1:TT
catch1 = catch
#catch1[4:TT] = rep(catch[4], TT-4)

#
plot(time, catch)

#
X1 = (sapply(time[-TT], FUN=pInd, t1=time))
X1 = cbind(1, X1)
sp1 = lm(catch~X1-1)
#


#
X2 = (sapply(time[-TT], FUN=pDiff, t1=time))
X2 = cbind(1, X2)
sp2 = lm(catch~X2-1)
#
lines(time, X2%*%sp2$coef)

#
X3 = (sapply(time[-c(TT, TT-1)], FUN=pCube, t1=time))
X3 = cbind(1, time, X3)
U = ncol(X3)
sp3 = lm(catch~X3-1)
#
L = length(sp3$coeff)
timePred = seq(1, TT, 0.1)
X3Pred = sapply(1:(length(sp3$coeff)-2), FUN=pCube, t1=timePred)
X3Pred = cbind(1, timePred, X3Pred)
lines(timePred, X3Pred%*%sp3$coeff) #X3Pred[,-(TT-1)]%*%sp3$coef[-(TT-1)])

#
X4 = bs(time, knots=time[-TT], intercept=T)
sp4 = lm(catch~X4)

##
#X2 = sapply(time[-TT], FUN=pDiff, t1=time)
#sp2 = lm(catch~X2)
#tCont = seq(1, TT, 0.1)
#X2Pred = sapply(time[-TT], FUN=pDiff, t1=tCont)
###
##plot(time, catch)
##lines(tCont, cbind(1, X2Pred)%*%sp2$coefficients)
#
##
##sp3 = lm(catch~X2^2)
##coeff = c(sp3$coefficients[1], sp3$coefficients[-1]/2)
##lines(tCont, cbind(1, X2Pred)%*%coeff)
#
##
#
#
###
##Xi  = sapply(time[-TT]-1, FUN=pInt, ti=time)
##Xim = sapply(time[-TT]-1, FUN=pInt, ti=time-1)
##X3  = Xi-Xim
##X3  = cbind(time-(time-1), X3[,-(TT-1)]) 
###X3  = X3[-1,-ncol(X3)]
##sp3 = lm(catch~X3-1)
##coef = sp3$coef#[-length(sp3$coef)] #c(sp3$coefficients[1], sp3$coefficients[-1])
#####coeff[1] = catch[1]
##f = function(t){
##	x = cbind(1, t(sapply(time[-TT], FUN=pDiff, t1=t)))	
##	x[-TT]%*%coef
##}
##f = Vectorize(f)
####
###lines(time, X3%*%coef, col='blue')
###lines(time, cbind(1, X2)[,-TT]%*%coef, col='red')
###points(time[-1], sapply(time[-1], FUN=function(x){integrate(f, x-1, x)$value}), col='red')
#
##
#X4 = matrix(0, nrow=TT, ncol=TT)
#for(i in time){
#	X4[i, 1] = 1 #time[i]-(time[i]-1)
#	if(i-1<1){ next }
#	for(j in 1:(i-1)){
#		X4[i, j+1] = ((time[i]^2)/2-time[i]*time[j]) - ((time[i-1]^2)/2-time[i-1]*time[j])
#	}
#}
##
#sp4 = lm(catch~X4-1)
#sp5 = glm(catch~X4-1, family=gaussian(link=log))
##sp6 = glm(catch~X4-1, family=gaussian(link=log), weights=1/log(catch)^2)
##sp5 = lm(catch~X4-1, wieght=1/catch^2)
##sp5 = rma(catch~X4-1, catch, method="FE")
##sp4 = gls(catch1~X4-1, weights=varPower())
#coef = sp4$coef
#
#png('spiky.png')
#plot(time, X4%*%coef, 'l', col='blue', ylim=c(0, 1500))
#lines(time, cbind(1, X2)%*%sp4$coefficients) 
#points(time, catch, pch=19)
#
##
#ff = function(t){
#        x = cbind(1, t(sapply(time[-TT], FUN=pDiff, t1=t)))
#	x%*%coef
#}
#ff = Vectorize(ff)
#points(time, sapply(time, FUN=function(x){integrate(ff, x-1, x)$value}), col='red')
#dev.off()
#
##
#sig = diag(1/catch^2)
#ISig = solve(sig)
#coefWLS = solve(t(X4)%*%ISig%*%X4) %*% t(X4)%*%ISig%*%catch
##sp5 = lm(catch~X4-1, wieght=1/catch^2)
#
##
#makeX = function(ts, knots){
#	#	
#	x = c()
#	for(t in ts){
#		tauStar = ceiling(t)-1 #t-1 #floor(t)
#		x = rbind( x, c( 1, (t^2/2-t*knots)*(t>knots) - (tauStar^2/2-tauStar*knots)*(tauStar>=knots)) )
#		#B_0: 1 -or- t-tauStar
#	}
#	return(as.data.frame(x))
#	#Tk = length(knots)
#	#x = matrix(0, nrow=1, ncol=Tk)	
#	#for(j in 1:Tk){
#	#	#
#	#	x[j] = (t^2/2-t*knots[j])*(t>knots[j]) - (tauStar^2/2-tauStar*knots[j])*(tauStar>=knots[j])
#	#}
#	#return(x)
#	#Ts = length(ts)
#	#X = matrix(0, nrow=Ts, ncol=length(knots))
#	#for(i in 1:Ts){
#	#	
#	#        X[i, 1] = 1 #time[i]-(time[i]-1)
#	#        if(i-1<1){ i=i+1; next }
#	#        for(j in 1:(i-1)){
#	#                X4[i, j+1] = ((time[i]^2)/2-time[i]*time[j]) - ((time[i-1]^2)/2-time[i-1]*time[j])
#	#        }
#	#	i = i+1
#	#}
#}
#
####
###knots = time
###t = 2.5
###tauStar = ceiling(t)-1
###x = (t^2/2-t*knots)*(t>knots) - (tauStar^2/2-tauStar*knots)*(tauStar>=knots)
###
##ts = seq(0, TT, 0.1)
##xs = makeX(ts, time[-TT])
##lines(ts, predict(sp4, xs))



























