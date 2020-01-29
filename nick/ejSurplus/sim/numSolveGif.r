rm(list=ls())

#
library(rstan)
library(foreach)
library(rootSolve)
library(doParallel)
registerDoParallel(cores=24)
#
source('funcs.r')
#
#set.seed(1)

#
#SOLVE
#

#
write.table(t(c('M', 'tail.x', 'tail.y', 'median.x', 'median.y', 'mean.x', 'mean.y')), file='bh_1800_b0_3000_fm_0.9_M_050_300Dat.csv', row.names=F, col.names=F, quote=F, sep=',')
write.table(t(c('delta', 'CStar', 'q', 's', 'alpha', 'beta', 'hStar', 'FStar', 'M')) , file='bh_1800_b0_3000_fm_0.9_M_050_300Par.csv', row.names=F, col.names=F, quote=F, sep=',')
#for(M in seq(0.01, 0.3, 0.002)){
#for(M in seq(0.04, 0.30, 0.01)){
#M=0.05:0.2, fix it
#F*/m=0.75, 1.5
outs = foreach( M=seq(0.04, 0.30, 0.005) )%dopar%{
#outs = foreach( M=seq(0.05, 0.2, 0.005) )%dopar%{
	##
	#dom = seq(0.05, 0.30, 0.005)
	#ind = which(dom==M)%%24
	#dir = sprintf('thread%d', ind)
	#suppressWarnings(dir.create(dir))
	#system(sprintf('cp numSolve.stan %s', dir))
	#setwd(dir)
	#
	#        ?         ~          *         *          ~       
	#       750/3000, 900/3000, 1200/3000, 1500/3000, 1800/3000
	#Bs/B0= 0.25    , 0.3     ,  0.4     , 0.5      , 0.6      
	#F*/M=0.9
	#
	#        ?         ~          *         *          *       
	#       750/3000, 900/3000, 1200/3000, 1500/3000, 1800/3000
	#Bs/B0= 0.25    , 0.3     ,  0.4     , 0.5      , 0.6      
	#F*/M=1.5
	#
	#        x         *          *         *          *          ~          x
	#       750/3000, 900/3000, 1200/3000, 1500/3000, 1800/3000, 2400/3000, 2700/3000
	#Bs/B0= 0.25    , 0.3     ,  0.4     , 0.5      , 0.6      , 0.8      , 0.9
	#F*/M=1.5
	#
	#        ?         ?          *         *          *          * 
	#       750/3000, 900/3000, 1200/3000, 1500/3000, 1800/3000, 2700/3000
	#Bs/B0= 0.25    , 0.3     ,  0.4     , 0.5      , 0.6      , 0.9
	#F*/M=5
	#
	Fs = M*0.9#1.5#5
	hs = 1-exp(-Fs)
	start = c(302, 0.12)
	names(start) = c('Cs', 'gamma')
	par = data.frame(
		kappa=0,
		wk=0.001,
		wI=0.001,
		Bs=1800,#2400,#2700,#1200,#1500,#900 ,#1000,#
		B0=3000,#3000,#3000,#3000,#3000,#3000,#2000,#
		M=M,
		hs=hs
	)
	#
	CG = multiroot(f, start, parms=par,
		maxiter=10^5,
		rtol=.Machine$double.eps,
		atol=.Machine$double.eps,
		ctol=.Machine$double.eps
	)$root
	Cs = CG[1]
	gamma = CG[2]	
	
	#
	#DATA
	#
	
	#
	q     = 1
	wk    = 0.001
	kappa = 0
	gaMod = -1
	delta = 1-exp(-M)
	#
	ab    = abGet(hs, Cs, delta, gamma, wk, kappa)
	alpha = ab$a
	beta  = ab$b
	#
	Y = 30
	C = c(seq(170, 200, 2), seq(198, 172, -2))-100 #rep(200, Y)#*1000/1500
	mod = model(q, alpha, beta, delta, kappa, gamma, Y, C)
	sd = 0.001
	cpue = rlnorm(Y, mod$logI, sd)
	#
	D = list(N=Y, cpue=cpue, C=C, gamma=gaMod, wk=wk, kappa=kappa, delta=delta)	
	write.table(t(c(delta, Cs, q, sd, alpha, beta, hs, Fs, M)), file='bh_1800_b0_3000_fm_0.9_M_050_300Par.csv', row.names=F, col.names=F, quote=F, sep=',', append=T)
	
	#
	#MODEL
	#
	
	#
	#interface w/ stan
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())
	out = stan( file = "numSolve.stan",
	        data  = D,
	        iter  = 10^4,
	        cores = 1,
		chains= 4,
		seed  = 1,
	        init  = function(){list(
	                q     = q,
	                hStar = hs,
	                CStar = Cs,
	                #delta = delta,
	                s     = sd
	        )}
	)
	#
	samples = extract(out, permuted=T)
	MM = length(samples[[1]])
	P = length(samples)-1
	
	##
	#m = stan_model( file="numSolve.stan" )
	#optOut = optimizing( m,
	#        data = D,
	#	seed = 1,
	#        init = function(){list(
	#                q     = q,
	#                hStar = hs,
	#                CStar = Cs,
	#                delta = delta,
	#                s     = sd
	#        )}
	#)
	
	#
	#LOOK
	#
	
	#
	postR = matrix(NA, nrow=MM, ncol=Y)
	postB = matrix(NA, nrow=MM, ncol=Y)
	postW = matrix(NA, nrow=MM, ncol=Y)
	post  = matrix(NA, nrow=MM, ncol=Y)
	pred  = matrix(NA, nrow=MM, ncol=Y)
	for(m in 1:MM){
	        #modOut = model(samples$q[m], samples$alpha[m], samples$beta[m], samples$delta[m], 0, -1, Y, C)
	        modOut = model(samples$q[m], samples$alpha[m], samples$beta[m], delta, 0, -1, Y, C)
		postR[m,] = modOut$R
	        postB[m,] = modOut$B
	        postW[m,] = modOut$W
	        post[m,] = exp(modOut$logI)
	        pred[m,] = rlnorm(Y, modOut$logI, samples$s[m])
	}
	#
	postUpper = matrix(NA, nrow=1, ncol=Y)
	postLower = matrix(NA, nrow=1, ncol=Y)
	#
	predUpper = matrix(NA, nrow=1, ncol=Y)
	predLower = matrix(NA, nrow=1, ncol=Y)
	for(t in 1:Y){
	        #
	        postUpper[t] = quantile(post[,t], 0.025)
	        postLower[t] = quantile(post[,t], 0.975)
	        #
	        predUpper[t] = quantile(pred[,t], 0.025)
	        predLower[t] = quantile(pred[,t], 0.975)
	}
	
	#
	#pairs(out)
	#
	
	##
	##pdf(sprintf('./pics/dataB0%d_M%d.pdf', par$B0, M*100))
	#year=1:Y
	#plot(year, cpue, ylim=c(min(pred), max(pred)))
	#polygon(c(year, rev(year)), c(predUpper, rev(predLower)), col='grey', border=NA)
	#polygon(c(year, rev(year)), c(postUpper, rev(postLower)), col='dimgrey', border=NA)
	#lines(year, colMeans(pred), lwd=5)
	#points(year, cpue)
	##dev.off()
	
	#
	#LOOK
	#
	
	#
	h0 = 0
	#
	s0 = sGet(h0, delta) #sGet(h0, samples$delta)
	W0 = WGet(s0, wk, kappa)
	r0 = rGet(h0, s0, W0, wk)
	R0 = RGet(samples$alpha, samples$beta, gaMod, r0)
	B0 = BGet(s0, W0, R0, wk)
	#
	sh = sGet(samples$hStar, delta) #sGet(samples$hStar, samples$delta)
	Wh = WGet(sh, wk, kappa)
	rh = rGet(samples$hStar, sh, Wh, wk)
	Rh = RGet(samples$alpha, samples$beta, gaMod, rh)
	Bh = BGet(sh, Wh, Rh, wk)
	##
	#s0Head = sGet(h0, optOut$par['delta'])
	#W0Head = WGet(s0Head, wk, kappa)
	#r0Head = rGet(h0, s0Head, W0Head, wk)
	#R0Head = RGet(optOut$par['alpha'], optOut$par['beta'], gaMod, r0Head)
	#B0Head = BGet(s0Head, W0Head, R0Head, wk)
	##
	#shHead = sGet(optOut$par['hStar'], optOut$par['delta'])
	#WhHead = WGet(shHead, wk, kappa)
	#rhHead = rGet(optOut$par['hStar'], shHead, WhHead, wk)
	#RhHead = RGet(optOut$par['alpha'], optOut$par['beta'], gaMod, rhHead)
	#BhHead = BGet(shHead, WhHead, RhHead, wk)
	##
	#FHead = -log(1-optOut$par['hStar'])
	#MHead = -log(1-optOut$par['delta'])
	##
	Fmsy = -log(1-samples$hStar)
	EM   = -log(1-delta)#-log(1-samples$delta)
	#
	headx = median(Fmsy/EM) #FHead/MHead   #
	heady = median(Bh/B0)   #BhHead/B0Head #
	tailx = Fs/M
	taily = par$Bs/par$B0
	#
	#pdf(sprintf('./pics/ratioB0%d_M%d.pdf', par$B0, M*100))
	png(sprintf('./pics/bh_1800_b0_3000_fm_0.9_M_050_300Gif%3d.png', round(M*1000, 0)))
	#dev.new()
	smoothScatter(Fmsy/EM, Bh/B0,
		xlab='F*/M',
		ylab=expression(B[h]/B[0]),
		xlim=c(0, 23), 
		ylim=c(0,1),
		main=sprintf("M=%1.2f", M)
	)#xlim=c(min(Fmsy/EM), max(c(tailx, headx, Fmsy/EM))), ylim=c(min(c(Bh/B0, taily, heady)), max(c(Bh/B0, taily, heady))))
	arrows(tailx, taily, x1=headx, y1=heady, col='red', lwd=3, length=0.15)
	points(tailx, taily, pch=19, col='red')
	dev.off()

	#
	#SAVE
	#	
	
	##
	#save.image(sprintf('./pics/bh_1800_b0_3000_M2_040_300RDat%3d.RData', round(M*1000, 0)))
	#
	write.table(t(c(M, tailx, taily, headx, heady, mean(Fmsy/EM), mean(Bh/B0))), file='bh_1800_b0_3000_fm_0.9_M_050_300Dat.csv', row.names=F, col.names=F, quote=F, sep=',', append=T)
	#
	return( out )
}
save.image(sprintf('./bh_1800_b0_3000_fm_0.9_M_050_300Dat.RData'))
