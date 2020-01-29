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
bRats = seq(0.2, 0.6, 0.1)
dRats = seq(0.5, 3, 0.5)
Ms = seq(0.05, 0.215, 0.007) 
#
outs = list()
for(bRat in bRats){ 
	#
	outs[[as.character(bRat)]]=list();
	for(dRat in dRats){
		#
	       	write.table(t(c('M', 'tail.x', 'tail.y', 'median.x', 'median.y', 'mean.x', 'mean.y')), file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30DatSD.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',')
	        write.table(t(c('delta', 'CStar', 'gamma', 'q', 'alpha', 'beta', 'hStar', 'FStar', 'M')) , file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30ParSD.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',')
		#
		#for(M in Ms){
		outs[[as.character(bRat)]][[as.character(dRat)]] = foreach( M=Ms )%dopar%{
			#
			Fs = M*dRat
                	hs = 1-exp(-Fs)
			#
			B0 = 3000
                	Bs = B0*bRat
			#
			start = c(302, 0.12)
        		names(start) = c('Cs', 'gamma')
			par = data.frame(
                	        kappa=0,
                	        wk=0.001,
                	        wI=0.001,
                	        Bs=Bs,
                	        B0=B0,
                	        M=M,
                	        hs=hs
                	)
			#
			CG = tryCatch({
				out = multiroot(f, start, parms=par,
                                	maxiter=10^5,
                                	rtol=.Machine$double.eps,
                                	atol=.Machine$double.eps,
                                	ctol=.Machine$double.eps
                        	)$root
				}, error=function(err){ return(NA) }
			)
			#
			if( any(is.na(CG)) ){
				write.table(t(c(NA, CG[1], CG[2], NA, NA, NA, hs, Fs, M)), file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30ParSD.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',', append=T)
				return(NA)
				#next 
			}	
			Cs = CG[1]
			gamma = CG[2]
			#return(CG)			

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
        		ab = tryCatch({
					out = abGet(hs, Cs, delta, gamma, wk, kappa)
				}, error=function(err){ return(NA) }
			)
			if( any(is.na(ab)) ){ 
			write.table(t(c(delta, Cs, gamma, q, ab$a,ab$b, hs, Fs, M)), file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30ParSD.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',', append=T)
				return(NA)
				#next 
			}
        		alpha = ab$a
        		beta  = ab$b
        		#
        		Y = 30
        		C = c(seq(170, 200, 2), seq(198, 172, -2))-100 		#rep(200, Y)#*1000/1500
			#
        		mod = model(q, alpha, beta, delta, kappa, gamma, Y, C) 
			#tryCatch({
			#		out = model(q, alpha, beta, delta, kappa, gamma, Y, C)
			#	}, error=function(err){ return(NA) }
			#)
			#if( any(is.na(mod$logI)) ){ 
			#	write.table(t(c(delta, Cs, gamma, q, alpha, beta, hs, Fs, M)), file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30Par.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',', append=T)
			#	#return(NA)
			#	next 
			#}
			#
        		sd = 0.1 #0.001
        		cpue = rlnorm(Y, mod$logI, sd)
        		#
        		D = list(N=Y, cpue=cpue, C=C, gamma=gaMod, wk=wk, kappa=kappa, delta=delta)
        		write.table(t(c(delta, Cs, gamma, q, alpha, beta, hs, Fs, M)), file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30ParSD.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',', append=T)
			
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
        		        chains= 2, #4
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
		
			#
        		#LOOK
        		#

        		#
        		h0 = 0
        		#
        		s0 = sGet(h0, delta) 					#sGet(h0, samples$delta)
        		W0 = WGet(s0, wk, kappa)
        		r0 = rGet(h0, s0, W0, wk)
        		R0 = RGet(samples$alpha, samples$beta, gaMod, r0)
        		B0 = BGet(s0, W0, R0, wk)
        		#
        		sh = sGet(samples$hStar, delta) 			#sGet(samples$hStar, samples$delta)
        		Wh = WGet(sh, wk, kappa)
        		rh = rGet(samples$hStar, sh, Wh, wk)
        		Rh = RGet(samples$alpha, samples$beta, gaMod, rh)
        		Bh = BGet(sh, Wh, Rh, wk)
			#
        		Fmsy = -log(1-samples$hStar)
        		EM   = M 						#-log(1-delta)#-log(1-samples$delta)
        		#
        		headx = median(Fmsy/EM) 				#FHead/MHead   #
        		heady = median(Bh/B0)   				#BhHead/B0Head #
        		tailx = Fs/M
        		taily = par$Bs/par$B0
        		#	
			png(sprintf('./pics/bRat_%1.2f_dRat_%1.1f_M_0.04_0.30Gif%1.3fSD.png', bRat, dRat, M))
        		#dev.new()
			plot(rep(dRats[1], length(bRats)), bRats, 
				xlim=c(0, 20), xaxs="i",
				ylim=c(0, 1),
				main=sprintf("M=%1.2f", M),
				xlab=expression(F[msy]/M),
                                ylab=expression(B[msy]/B[0])
			)
			for(d in dRats[-1]){ points(rep(d, length(bRats)), bRats) }
        		smoothScatter(Fmsy/EM, Bh/B0,
                	        xlab=expression(F[msy]/M),
                	        ylab=expression(B[msy]/B[0]),
                	        xlim=c(0, 20), xaxs="i",
                	        ylim=c(0,0.5),
                	        main=sprintf("M=%1.2f", M),
				nrpoints=0,
				add=T
                	)	
        		arrows(tailx, taily, x1=headx, y1=heady, col='red', lwd=3, length=0.15)
        		points(tailx, taily, pch=19, col='red')
        		dev.off()
			
                	#
                	#SAVE
                	#	
                	
                	#
                	write.table(t(c(M, tailx, taily, headx, heady, mean(Fmsy/EM), mean(Bh/B0))), file=sprintf('bRat_%1.2f_dRat_%1.1f_M_0.04_0.30DatSD.csv', bRat, dRat), row.names=F, col.names=F, quote=F, sep=',', append=T)
			#
			return( out )
		}
		#
        	save.image(sprintf('./bRat_%1.2f_dRat_%1.1f_M_0.04_0.30DatSD.RData', bRat, dRat))
        	system(sprintf('convert -delay 20 ./pics/bRat_%1.2f_dRat_%1.1f_M_0.04_0.30Gif*SD.png ./pics/bRat_%1.2f_dRat_%1.1f_M_0.04_0.30SD.gif', bRat, dRat, bRat, dRat))
	}
}
