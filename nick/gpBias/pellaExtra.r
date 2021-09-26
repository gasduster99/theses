rm(list=ls())

#
library(VGAM)
library(boot)
library(pracma)
library(mvtnorm)
library(graphics)
library(parallel)
library(rootSolve)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)
#
source('gpFunk.r')
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
        #linearly interpolate catches
	ft = floor(t)
	q  = (t-ft) 
        Cl = catch[ft]
	Cu = catch[ft+1]
	C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
	#
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
        out = R - C
        #
        return( list(out) )
}

#
getPar = function(xi, zeta, M){
        #
        gamma = 1/zeta
        odd = (1-zeta)/zeta
        alpha = xi*M*odd*(1-zeta)^(-odd)
        #
        return(c(alpha, gamma))
}

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
#P0Funk = function(alpha, beta, gamma, M){ (alpha/(M)-1)^gamma * beta^-gamma }
#exp(gamma*log(alpha/(M+ff)-1) - gamma*log(beta)) }
#PBar = function(ff, alpha, beta, gamma, M){ (alpha/(M+ff)-1)^gamma * beta^-gamma }
PBar = function(ff, alpha, beta, gamma, M){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }

##
#getBeta = function(alpha, gamma, M, P0){
#        optim(0.01, function(x){ (PBar(0, alpha, x, gamma, M)-P0)^2 }, method='Brent', lower=eps(), upper=10^6)$par
#}
#
##
#FMsy = function(alpha, gamma, M){
#        #
#        FUpper = alpha-M
#        root = uniroot.all(function(ff){ 1 - exp(log(ff)+log(alpha)+log(gamma) - (2*log(M+ff)+log(alpha/(M+ff)-1))) }, c(0, FUpper))
#        #
#        if(length(root)<1){ return(NA) }
#        if(length(root)>1){ return(root[1]) }
#        #
#        return(root)
#}

#
FMsy = function(alpha, gamma, M){
        ##
        #FUpper = P0 #alpha-M
        ##root = uniroot.all(function(ff){  ((gamma-1)/alpha)^(gamma-1) - (M+ff)^(1/(gamma-1)) - (ff/(
        #root = uniroot.all(function(ff){ 1 - ((M+ff)*((gamma-1)/alpha))^(1/(gamma-1)) - (ff/(gamma-1)
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        #
        root = alpha/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
        return(root)
}


##
#dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/(x+2))^2) }
dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-0.5)^2) }

###
##getData = function(dir, xiSims, zetaSims){
##	#dir	: a directory containing data
##	#xiSims	: xis of simulated data
##	#zetaSims: zetas of simulated data
##	
##	#
##	di = 1
##	D = data.frame(xi=double(), zeta=double(), xiHat=double(), xiMin=double(), zetaHat=double(), distHat=double(), distMin=double(), cllhs=double(), cllhsV=double(), stringsAsFactors=F)
##	for(i in 1:length(zetaSims)){
##        	for(j in 1:length(xiSims)){
##                	#
##                	fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
##                	fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
##			
##			#
##			if( file.exists(fileDat) & file.exists(fileFit) ){
##                		#
##                		dat = readRDS(fileDat)
##                		fit = readRDS(fileFit)
##
##               		#
##                		xiHat   = -log(1-cloglog(fit$cllhs, inverse=T))/M
##                		zetaHat = 1/(xiHat+2)
##                		distHat = sqrt((dat$xi-xiHat)^2 + (dat$zeta-zetaHat)^2)
##                		opt = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)
##                		#
##				if(length(fit$rsCov)==0){ v=0 }else{ v=fit$rsCov['cllhs', 'cllhs'] }
##                		D[di,]  = c(dat$xi, dat$zeta, xiHat, opt$minimum, zetaHat, distHat, opt$objective, fit$cllhs, v)
##                		di=di+1
##                	}
##		}
##	}
##	#
##	return(D)
##}
#
##
#getData = function(dir, xiSims, zetaSims){
#	#dir	: a directory containing data
#	#xiSims	: xis of simulated data
#	#zetaSims: zetas of simulated data
#	#
#	#Seed values are the initiating values to launch the inversion
#	#Inv values are the numerical inversions actually found for the 3 parameter curve at each seed value
#	#BH values are the MLE fits under the BH model.
#	
#	#
#	di = 1
#	D = data.frame(xiSeed=double(), zetaSeed=double(), xiInv=double(), zetaInv=double(), xiBH=double(), zetaBH=double(), minDist=double(), lF=double(), lFV=double(), stringsAsFactors=F)
#	for(i in 1:length(zetaSims)){
#        	for(j in 1:length(xiSims)){
#                	#
#                	fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
#                	fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
#			
#			#
#			if( file.exists(fileDat) & file.exists(fileFit) ){
#				#
#				dat = readRDS(fileDat)
#				fit = readRDS(fileFit)	
#				
#				#the inversion actually found
#				FInv = FMsy(dat$alpha, 1, M)
#				xiInv = FInv/M
#				zetaInv = PBar(FInv, dat$alpha, dat$beta, dat$gamma, M)/PBar(0, dat$alpha, dat$beta, dat$gamma, M)
#				#the bh fit found
#				Fs = FMsy(fit$alpha, 1, M)
#				xiBH   = Fs/M 
#				zetaBH = PBar(Fs, fit$alpha, fit$beta, fit$gamma, M)/PBar(0, fit$alpha, fit$beta, fit$gamma, M)
#				md = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)$objective				
#
#				#NOTE: replace with a propper observation uncertainty
#				if(length(fit$rsCov)==0){ v=0 }else{ v=getlFV(fit)$lFV }
#				#print( c(dat$xi, dat$zeta, xiHat, zetaHat, md, log(Fs), v) )
#				D[di,] = c(dat$xi, dat$zeta, xiInv, zetaInv, xiBH, zetaBH, md, log(Fs), v)
#				#print(dim(D))
#				di = di+1
#			}
#		}
#	}
#	#
#	return(D)
#}

##
#getlFV = function(fit, MM=10^4, samples=F){
#	#
#	who = c('lalpha')
#	C = fit$rsCov[who, who]
#	m = c(fit$lalpha)
#	sam = rnorm(MM, m, sqrt(C)) #rmvnorm(M, m, C)
#	#
#	als = exp(sam)
#	lFs = sapply(als, function(a){log(FMsy(a, 1, M))})
#	#
#	out = list(lF=mean(lFs, na.rm=T), lFV=var(lFs, na.rm=T))
#	if( samples ){ out$lFSamples=lFs }
#	#
#	return(out)
#}

#
#
#

#
#dat = readRDS('modsPellaFineQFixReduxP06000/datGen_xi2.95_zeta0.53.rda')
#fit = readRDS('modsPellaFineQFixReduxP06000/fit_xi2.95_zeta0.53.rda')
##
#dat = readRDS('modsPellaFineQFixReduxP06000/datGen_xi1_zeta0.15.rda')
#fit = readRDS('modsPellaFineQFixReduxP06000/fit_xi1_zeta0.15.rda')

#
M = 0.2
P0 = 10000
catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 
254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
TT = length(catch)

#
xi = 3.5
zeta = 0.2
ag = getPar(xi, zeta, M)
alpha = ag[1]
gamma = ag[2]
beta = P0
##
#datGen = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #model
#        time=c(1:TT), catch=catch[1:TT], M=M,                    #constants
#        alpha=alpha, beta=beta, gamma=gamma,            #parameters
#        lalpha=log(alpha), lbeta=log(beta),             #reparameterize
#        lq=log(0.00049), lsdo=log(0.01160256) #log(0.1160256)           #nuisance paramete
#        #xi=xi, zeta=zeta                                #other incidentals to carry along
#)
#datGen$iterate("lsode") #"vode")#("euler")#("daspk")
#
##
#cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
#
##
#fit = prodModel$new(
#        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #model
#        time=1:TT, catch=catch, M=M,                            #constants
#        alpha=alpha, beta=P0, gamma=2,  #parameters
#        lalpha=log(alpha), lbeta=log(P0),       #reparameterize
#        lq=log(0.00049), lsdo=log(0.01160256)#, #log(0.1160256),                 #nuisance parameters
#        #xi=xiSims[j], zeta=zetaSims[i]                          #other incidentals to carry along
#)
#fit$iterate("vode")

##
#optAns = fit$optimize(cpue,
#        c('lsdo', 'lalpha'), #'lq'),
#        lower   = c(log(0.001), log(M)), #log(1e-7)),
#        upper   = c(log(1), log(100)), #log(1e-2)),
#        gaBoost = list(run=10, parallel=FALSE, popSize=10^3),
#        persistFor=5
#)

#
dir = "./modsPellaFineQFixReduxP010000"
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
cpue = rlnorm(TT, dat$lq+log(dat$N), exp(dat$lsdo))
#
fit = readRDS(fileFit)
fit$iterate("vode")

#
optAns = fit$optimize(cpue,
        c('lsdo', 'lalpha', 'lbeta'),
        lower   = c(log(0.001), log(0), log(P0/2)),     #log(10^3)), 
        upper   = c(log(1), log(100), log(2*P0)),       #log(10^4)),
        gaBoost = list(run=100, parallel=6, popSize=10^3),
        persistFor=0, cov=T #persistFor=5, cov=T
)
#optAns = fit$optimize(cpue,
#        c('lsdo', 'lalpha', 'lbeta'),
#        lower   = c(log(0.001), log(M), log(P0/2)),
#        upper   = c(log(1), log(100), log(2*P0)),
#        cov     = T
#        #c('lsdo', 'alpha', 'beta', 'lq'),                                                      
#        #lower   = c(log(0.001), eps(), eps(), log(1e-7)),
#        #upper   = c(log(1), 10, 2, log(1e-2)),
#)






##
##DATA
##
#
##
#dir = "./modsShepFineQFix/" #"./modsShepFine/" #"./modsShepTry/" #'./modsFine/' #'./modsHess/'
##
#zetaSims = seq(0.1, 0.8, 0.01) #seq(0.1, 0.8, 0.05) 	#rev(seq(0.1, 0.80, 0.01)) #
#xiSims =   seq(0.5, 3.5, 0.05) #seq(0.5, 3.5, 0.25)		#seq(0.5, 3.5, 0.05)       #
#
##
#M = 0.2
##time: 70
#D = getData(dir, xiSims, zetaSims)
##
#bub = 0.2
#D = D[D$xiInv<=max(xiSims)*(1+bub),]
#D = D[D$xiInv>=min(xiSims)*(1-bub),]
#D = D[D$zetaInv<=max(zetaSims)*(1+bub),]
#D = D[D$zetaInv>=min(zetaSims)*(1-bub),]
#D = D[D$lFV!=0 & D$xiBH<20,] #lalpha==0.04280697 is a numerical issue
##
#png('shepDat.png')
#plot(D[,c("xiInv", "zetaInv")], 
#	ylim=c(min(D[,c("zetaInv", "zetaBH")]), max(D[,c("zetaInv", "zetaBH")])), 
#	xlim=c(min(D[,c("xiInv", "xiBH")]), max(D[, c("xiInv", "xiBH")])),
#	main='Shepherd',
#	col=map2color(D$minDist, hcl.colors(60, "Zissou 1", rev=T)),
#	pch=20
#)
#points(D[,c("xiBH", "zetaBH")], col=map2color(D$minDist, hcl.colors(60, "Zissou 1", rev=T)))
#dev.off()
#
##
##
##
#
##pick a polynomial mean function
#Tg = diag(D$lFV)
#X = cbind(1, D$xiInv, D$zetaInv)
#axes = X[,2:3]
#registerData(D$lF, X, axes, Tg)
#par = c(l1=0.5, l2=0.5, th=eps(), nu=1, s2=0.1)
#gpFit = gpMAP(par, hessian=F, psiSample=F)
#print(gpFit)
#
##prediction
#zetaStar = seq(min(D$zetaInv), max(D$zetaInv), 0.001) 	#seq(0.15, 0.35, 0.001)  #rev(seq(0.1, 0.80, 0.01)) #
#xiStar   = seq(min(D$xiInv), max(D$xiInv), 0.005)  	#seq(1, 3.5, 0.005)
#XStar = cbind(1, expand.grid(xiStar, zetaStar))
#gpPred = gpPredict(XStar, XStar[,2:3], gpFit)
#
##bias
#xiHat = exp(gpPred)/M
#xBias = sweep(xiHat, 1, xiStar)
#yBias = sweep(1/(xiHat+2), 2, zetaStar)
##
#eucBias = mcmapply(function(xiHat, xi, zeta){
#                dist(xiHat, xi, zeta)
#        }, xiHat, XStar[,2], XStar[,3], mc.cores=detectCores()
#)
#eucBias = matrix(eucBias, nrow=length(xiStar), ncol=length(zetaStar))
#
##xi bias
#png("xBias.png")
##
#maxXBias = abs(max(xBias, na.rm=T))
#minXBias = abs(min(xBias, na.rm=T))
#posCols = hcl.colors(round(100*maxXBias/(maxXBias+minXBias)), "Reds 2", rev=T)
#negCols = hcl.colors(round(100*minXBias/(maxXBias+minXBias)), "Blues 2", rev=F)
#xCols = c(negCols, "#FFFFFF", posCols)
##
#par(mar=c(5, 4, 4, 5)+0.1)
#image(xiStar, zetaStar, xBias,
#	col  = adjustcolor(xCols, alpha.f=0.6),  #hcl.colors(41, "RdBu", rev=T)
#        xlab = 'Xi',
#        ylab = 'Zeta',
#	main = "Bias in Estimated Optimal Fishing"
#)
#curve(1/(x+2), from=0, to=4, lwd=3, add=T) #col=map2color(0, hcl.colors(41, "RdBu", rev=T)),
#show = seq(1, length(xCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.1f", rev(seq(min(xBias, na.rm=T), max(xBias, na.rm=T), length.out=length(show)))),
#        fill = rev(xCols[show]), #colMap[c(1, 10, 20)], 
#        xpd = NA
#)
#dev.off()
#
##zeta bias
#png("yBias.png")
##
#maxYBias = abs(max(yBias, na.rm=T))
#minYBias = abs(min(yBias, na.rm=T))
#posCols = hcl.colors(round(100*maxYBias/(maxYBias+minYBias)), "Reds 2", rev=T)
#negCols = hcl.colors(round(100*minYBias/(maxYBias+minYBias)), "Blues 2", rev=F)
#yCols = c(negCols, "#FFFFFF", posCols)
##
#par(mar=c(5, 4, 4, 5)+0.1)
#image(xiStar, zetaStar, yBias,
#	col  = adjustcolor(yCols, alpha.f=0.6), 
#        xlab = 'Xi',
#        ylab = 'Zeta',
#	main = "Bias in Estimated Optimal Biomass"
#)
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#show = seq(1, length(yCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.2f", rev(seq(min(yBias, na.rm=T), max(yBias, na.rm=T), length.out=length(show)))),
#        fill = rev(yCols[show]), #colMap[c(1, 10, 20)], 
#        xpd = NA
#) 
#dev.off()
#
##euc bias
#png("directionalBias.png")
##
#eucCols = hcl.colors(41, "Reds 2", rev=T)
##
#par(mar=c(5, 4, 4, 5)+0.1)
#image(xiStar, zetaStar, eucBias,
#        col  = adjustcolor(eucCols, alpha.f=0.6),
#        xlab = 'Xi',
#        ylab = 'Zeta',
#	main = "Directional Bias"
#)
#curve(1/(x+2), from=0, to=4, lwd=3, add=T) 
#w = (XStar[,2]>0.5 & XStar[,2]<3.5 & XStar[,3]>0.2 & XStar[,3]<0.75) 
#thin = c(T,rep(F,60))
#quiver(
#        XStar[w,2][thin], XStar[w,3][thin],
#        xBias[w][thin], yBias[w][thin],
#        scale=0.025
#)
#show = seq(1, length(eucCols), length.out=20)
#legend(grconvertX(415, "device"), grconvertY(90, "device"), #grconvertX(0.5, "device"), grconvertY(1, "device"),  #
#        sprintf("%1.1f", rev(seq(min(eucBias, na.rm=T), max(eucBias, na.rm=T), length.out=length(show)))),
#        fill = rev(eucCols[show]), #colMap[c(1, 10, 20)], 
#        xpd = NA
#)
#dev.off()
#
##residuals






#
#OLD
#

#X = cbind(D$xi, D$zeta)
#
##
#xiPred   = seq(0.7,3.5,length.out=50)	#seq(0.5, 3.5, 0.01) 		#	
#zetaPred = seq(0.2,0.75,length.out=50)	#seq(0.1, 0.75, 0.01)   	#
#predMesh = as.matrix(expand.grid(xiPred, zetaPred))
#colnames(predMesh) = c('xi', 'zeta')
##
#line = lm(cllhs~xi+zeta, data=D)
##gp = gpModel$new( S2=S2, X=X, Y=D$cllhs, obsV=D$cllhsV, #lm=line, 
##        B  = line$coeff,	#c(4.1894, -0.1452, -10.3066), 
##	v0 = 0.343979394716728,
##        v1 = 0.0270240223433852,
##        s2 = 2.93501281,
##        cr = -0.433574770896697
##)
###
##optOut0 = gp$fit(c('v0', 'v1', 's2', 'cr'),
##        lower   = c(eps(), eps(), eps(), -1),
##        upper   = c(Inf, Inf, Inf, 1),
##        cov     = T
##)
##gp$save('gpHotStart.rda')
#gp = readRDS('gpHotStart.rda')
#gp$printSelf()
##
#bnd = rep(Inf, length(line$coeff))
#optOut1 = gp$fit(c('B'),
#        lower   = c(-bnd),
#        upper   = c(bnd),
#        cov     = T
#)
#gp$save('gpHotterStart.rda')
#gp$printSelf()
##
#bnd = rep(Inf, length(line$coeff))
#optOut2 = gp$fit(c('v0', 'v1', 's2', 'cr', 'B'),
#        lower   = c(eps(), eps(), eps(), -1, -bnd),
#        upper   = c(Inf, Inf, Inf, 1, bnd),
#	#gaBoost = list(run=30, popSize=detectCores(), parallel=detectCores()),
#        cov     = T
#)
#gp$save('gpBigOpt.rda')
#gp$printSelf()
#
##
#yy = gp$predictMean(predMesh)
##SS = diag(gp$predictVar(predMesh))
###
##xS = exp(2*yy+SS)*(exp(SS)-1)/M^2
##K = 10^6
##yS = mcmapply(function(yy, SS){
##		logitM(2, log(2*M)-yy, SS, K=K) - logitM(1, log(2*M)-yy, SS, K=K)^2
##	}, yy, SS, mc.cores=8
##)
##xS = matrix(xS, nrow=length(xiPred))
##yS = matrix(yS, nrow=length(xiPred))
#
##
#xDist = matrix((exp(yy)/M)-xiPred, nrow=length(xiPred)) #, byrow=T)
#yDist = matrix((1/(exp(yy)/M+2))-zetaPred, nrow=length(xiPred))#, byrow=T)
##
#bigD = mcmapply(function(xiHat, xi, zeta){
#                dist(xiHat, xi, zeta)
#        }, exp(yy)/M, xiPred, zetaPred, mc.cores=detectCores()
#)
#
##
##OUTPUT
##
#
##
#m = 2.5
#
##
#pdf('gpXiBiasFineRB.pdf')
##dev.new()
#filled.contour(xiPred, zetaPred, xDist, 
#	zlim=c(-m, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#		#
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#	},
#	plot.title = title(
#		main = "Bias in Estimated Optimal Fishing Rate",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#pdf('gpZetaBiasFineRB.pdf')
##dev.new()
#filled.contour(xiPred, zetaPred, yDist, 
#	zlim=c(-0.55,0.55),
#        ylab = expression(frac(F[msy]/F[0])),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')	
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#	},
#	plot.title = title(
#		main = "Bias in Estimated Optimal Biomass",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3"),
#)
#dev.off()
#
##
#pdf('gpBiasArrowB.pdf')
##dev.new()
#datGrid = expand.grid(xiSims, zetaSims)
#colnames(datGrid) = c('xi', 'zeta')
#w = (predMesh[,'xi']>0.5 & predMesh[,'xi']<3.5 & predMesh[,'zeta']>0.2 & predMesh[,'zeta']<0.75) #paste0(predMesh[,'xi'], predMesh[,'zeta'])%in%paste0(datGrid[,'xi'], datGrid[,'zeta']) #
#thin = c(T,rep(F,3))
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred)),
#        zlim=c(0, m),
#        plot.axes = {
#                quiver(
#			predMesh[w,'xi'][thin], predMesh[w,'zeta'][thin],
#			xDist[w][thin], yDist[w][thin],
#			scale=0.05
#		) 
#                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#
#		axis(1, seq(min(xiPred), max(xiPred), by=0.5))
#                axis(2, seq(min(zetaPred), max(zetaPred), by=0.1))
#        },
#	plot.title = title(
#		main = "Directional Bias",
#		ylab = expression(B[msy]/B[0]),
#        	xlab = expression(F[msy]/M),	
#	),
#	#key.title = title(main = "Directional Bias"),
#        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
#)
#dev.off()


##
#source("gpClass0.0.1.r")



##X's are 2-vectors; L is 2x2; s2 is scalar
#S2 = function(X0, X1, s2, v0, v1, cr){
#	maxD = dbinorm(X0[,1], X0[,2], X0[,1], X0[,2], v0, v1, sqrt(v0*v1)*cr, log=T)
#	s2*mcmapply(function(x01, x02, m){
#		exp(dbinorm(X1[,1], X1[,2], x01, x02, v0, v1, sqrt(v0*v1)*cr, log=T)-m)
#	}, X0[,1], X0[,2], maxD, mc.cores=detectCores())	
#}
#
##
#logitM = function(n, mu, sig, K=10^5){
#	mean( inv.logit(qnorm((1:(K-1))/K, mu, sig))^n )
#}













#
#dev.new()
#filled.contour(xiPred, zetaPred, yDist/sqrt(yS), 
#	zlim=c(-0.6,0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)




##all parameters scalar
#Kr = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
#K = function(x0, x1, l){ exp((-(x0-x1)^2)/(2*l^2)) }
#Kmat = function(X, l){ 
#	#X	:a vector of predictors
#	#l	:a scalar length scale
#	
#	#
#	sapply(X, function(x){K(x,X,l)})
#}
##
#Kn = function(X, Y, lx, ly){ 
#	#X	:a vector of predictors
#	#Y	:a vector of predictors
#	#lx	:a scalar length scale for X
#	#ly	:a scalar length scale for Y
#	
#	#
#	sapply(Y, function(y){K(y,Y,ly)})*sapply(X, function(x){K(x,X,lx)})
#}
##
#Ku = function(X, Y, lx, ly){ 
#	#X	:a vector of predictors
#	#Y	:a vector of predictors
#	#lx	:a scalar length scale for X
#	#ly	:a scalar length scale for Y
#	
#	#
#	sapply(Y, function(y){K(y,Y,ly)})+sapply(X, function(x){K(x,X,lx)})
#}
##
#KnPred = function(Xs, Ys, X, Y, lx, ly){
#	sapply(Y, function(y){K(y,Ys,ly)})*sapply(X, function(x){K(x,Xs,lx)})
#}
#
##
#loglikeL = function(y, D, lx, ly, r, s2){	#function(y, D, L, sig2){ #
#	#
#	X = cbind(D$xi, D$zeta)
#	SIG2 = S2(X, X, lx, ly, r, s2) + diag(D$cllhsV)
#		
#	#
#	return( dmvnorm(y, sigma=SIG2, log=T) )
#}





#pp = predict( line )
#ppPred = predict(line, predMesh)


###
##optOut = c(mean(c(0.34395317, 0.02705277)), 2.93501281) 
##names(optOut) = c('lx', 'sig2')
#optOut = c(0.34395317, 0.02705277, -0.43195056, 2.93501281) 
#names(optOut) = c('lx', 'ly', 'cr', 'sig2')
##optOut = c(0.34395317, 0.02705277, 2.93501281)
##names(optOut) = c('lx', 'ly', 'sig2') 
##time: 340
#t = system.time({ 
#	optOut = optim(optOut, 
#		function(x){ 	
#			-loglikeL(D$cllhs-pp, D, x[1], x[2], x[3], x[4])
#		},
#		lower = c(10^-6, 10^-6, -1, 0),
#		upper = c(Inf, Inf, 1, Inf),
#		method = "L-BFGS-B",
#		hessian = T
#	)
#})
#writeLines(sprintf("Opt: %s",t[3]))
##IL = matrix(c(optOut$par[1], optOut$par[3], optOut$par[3], optOut$par[2]), nr=2)
#lx   = optOut$par['lx']
#ly   = optOut$par['ly']
#r    = optOut$par['cr']
#sig2 = optOut$par['sig2']
## 
##KInv = chol2inv(chol( sig2*Kn(D$xi,D$zeta,lx,ly) + diag(D$cllhsV) ))
#X = cbind(D$xi, D$zeta)
#KInv = chol2inv(chol( S2(X, X, lx, ly, r, sig2) + diag(D$cllhsV) ))
##KPX  = sig2*KnPred(predMesh$xi, predMesh$zeta, D$xi, D$zeta, lx, ly)
#KPX = S2(X, predMesh, lx, ly, r, sig2)
##KXP  = sig2*KnPred(D$xi, D$zeta, predMesh$xi, predMesh$zeta, lx, ly)
##KPP  = sig2*KnPred(predMesh$xi, predMesh$zeta, predMesh$xi, predMesh$zeta, lx, ly)
####time: 66
#ys = ppPred + KPX%*%KInv%*%(D$cllhs-pp) #cloglog(KPX%*%KInv%*%(D$cllhs-pp), inverse=T) #
###time: 470
##Ss = KPP-KPX%*%KInv%*%KXP
###writeLines(sprintf("Krig Var: %s",t[3]))

##
##CLLHS OUTPUT
##
#
##
#pdf('gpCllhsFineR.pdf')
#filled.contour(xiPred, zetaPred, matrix(ys, nrow=length(xiPred), byrow=T), 
#	zlim=c(-3.5, 1.5), #c(-3,0),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	}
#)
#dev.off()
#
##
#m = 3.5
#
##
#pdf('gpXiBiasFineR.pdf')
#xDist = matrix((exp(ys)/M)-xiPred, nrow=length(xiPred), byrow=T)
##xDist[xDist<(-1.2)] = -10
#filled.contour(xiPred, zetaPred, xDist, 
#	zlim=c(-m, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#pdf('gpZetaBiasFineR.pdf')
#yDist = matrix((1/(exp(ys)/M+2))-zetaPred, nrow=length(xiPred), byrow=T)
#filled.contour(xiPred, zetaPred, yDist, 
#	zlim=c(-0.6,0.6),
#	plot.axes = {
#		points(D$xi, D$zeta)
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blue-Red 3")
#)
#dev.off()
#
##
#bigD = mapply(function(xiHat, xi, zeta){ 
#		dist(xiHat, xi, zeta) 
#	}, exp(ys)/M, xiPred, zetaPred
#) 
##
#pdf('gpBiasFineR.pdf')
##m = 3.5 #quantile(abs(bigD), 0.5)
#filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T), 
#	zlim=c(0, m),
#	plot.axes = {
#		points(D$xi, D$zeta, pch='.')
#		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
#		#points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
#		#points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
#	},
#	color.palette = function(n) hcl.colors(n, "Blues", rev = TRUE)
#)
#dev.off()


###
##w = predMesh$X>0.5 & predMesh$X<3.5 & predMesh$Y>0.2 & predMesh$Y<0.75
###
##pdf('gpBiasArrow.pdf')
##filled.contour(xiPred, zetaPred, matrix(bigD, nrow=length(xiPred), byrow=T),
##        zlim=c(0, m),
##        plot.axes = {
##                quiver(predMesh$X[w[c(T,F)]], predMesh$Y[w[c(T,F)]], xDist[w[c(T,F)]], yDist[w[c(T,F)]])
##                #for(i in 1:length(predMesh$xi)){
##                #       arrows(predMesh$xi[i], predMesh$zeta[i], predMesh$xi[i] + scale * u, predMesh$zeta 
##                #}
##                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
##        },
##        color.palette = function(n) hcl.colors(n, "Reds", rev = TRUE)
##)
##dev.off()
#
#
#
#
#
##sapply(1:length(xiPred),
##	function(xi){sapply(1:length(zetaPred), 
##		function(zeta){ dist((exp(ys)/M)[], xiPred[], zetaPred[]) }
##	)
##}
#
###
##pdf('gpCllhsHessVar.pdf')
##filled.contour(xiPred, zetaPred, matrix(diag(Ss), nrow=length(xiPred), byrow=T), 
##	zlim=c(0, 0.6),
##	plot.axes = {
##		points(D$xi, D$zeta)
##		lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3)
##	}
##)
##dev.off()
#
#
#
#
###
##p = scatterplot3d(D$xi, D$zeta, D$cllhs)
##plane = lm(cllhs~xi+zeta, data=D)
##p$plane3d(plane)
##pp = predict(plane)
##scatterplot3d(D$xi, D$zeta, D$cllhs-pp)
#
###
###DISTANCE OUTPUT
###
##
###
##plane = lm(log(distHat-distMin)~xi+zeta, data=D)
##pp = predict(plane)
##yy = log(D$distHat-D$distMin)-pp
###
##yys = KPX%*%KInv%*%(yy)
###
##pdf('gpLogBiasNoHess.pdf')
##filled.contour(xiPred, zetaPred, matrix(yys, nrow=length(xiPred), byrow=T),
##        #zlim=c(-50, 20),
##        plot.axes = {
##                points(D$xi, D$zeta)
##                lines(seq(0,4,0.1), 1/(seq(0,4,0.1)+2), lwd=3) 
##                #points(D$xi[D$cllhsV==0], D$zeta[D$cllhsV==0], pch=16)
##                #points(D$xi[D$cllhsV==max(D$cllhsV)], D$zeta[D$cllhsV==max(D$cllhsV)], pch=16, col='blue')
##        }
##)
##dev.off()
