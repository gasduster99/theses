rm(list=ls())

#
library(rstan)
library(foreach)
library(doParallel)
library(RColorBrewer)

#
#FUNCTIONS
#

#
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

#
#DATA
#

#
MS = seq(0.06, 0.30, 0.005)
glob = Sys.glob('dRat_*_M_0.04_0.30Dat.RData')
dRats = c()
for(g in glob){ dRats = c(dRats, as.numeric(strsplit(g, "_")[[1]][2])) }
bRats = dRats/(dRats*(dRats+2))
#
forOut = list() #bList=list(), dList=list())
#registerDoParallel(cores=24)
#forOut = foreach( g=glob )%dopar%{
for(g in glob){
	#
	forOut[[g]] = list()
	forOut[[g]][['dList']] = list()
        forOut[[g]][['bList']] = list()
	##
	#forOutie = list()
	#forOutie[['dList']] = list()
        #forOutie[['bList']] = list()
	##
	load(g)
	O = length(outs)
	for(o in 1:O){
		#
		M = MS[o]	
		out = outs[o][[1]]
		#
		wk    = 0.001
		delta = 1-exp(-M)
		kappa = 0
		gaMod = -1
		#
		samples = extract(out, permuted=T)
        	MM = length(samples[[1]])
        	P = length(samples)-1
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
		#
		Fmsy = -log(1-samples$hStar)
                EM   = -log(1-delta) #-log(1-samples$delta)
		#
		forOut[[g]][['dList']][[as.character(M)]] = Fmsy/EM
		forOut[[g]][['bList']][[as.character(M)]] = Bh/B0
                ##
		#forOutie[['dList']][[as.character(M)]] = Fmsy/EM
		#forOutie[['bList']][[as.character(M)]] = Bh/B0
       	}
	##
	#return( forOutie )
}

#
#PLOT
#

#
cols = rep(brewer.pal(11, 'Spectral')[-6], 2)
EM = names(forOut[[g]][['dList']])
#
for(em in EM){
	png(sprintf('./pics/M_0.04_0.30CloudGif%1.3f.png', as.numeric(em)))
	smoothScatter(forOut[[1]][['dList']][[em]], forOut[[1]][['bList']][[em]],
		xlab=expression(F[msy]/M),
		ylab=expression(B[msy]/B[0]),
		xlim=c(0, 3.5), xaxs="i",
		ylim=c(min(bRats)-0.05,0.5), #yaxs="i",
		main=sprintf("M=%1.2f", as.numeric(em)),
		colramp=colorRampPalette(c(addalpha('white', 0), addalpha(cols[1], 0.9)), alpha=T),
		nrpoints=0
	)
	points(dRats[1], bRats[1], col=cols[1], pch=19)
	#
	i = 2
	for(g in glob[-1]){
		#
		smoothScatter(forOut[[g]][['dList']][[em]], forOut[[g]][['bList']][[em]],
	        	xlab=expression(F[msy]/M),
	        	ylab=expression(B[msy]/B[0]),
	        	xlim=c(0, 3.5), xaxs="i",
	        	ylim=c(min(bRats)-0.1,0.5), #yaxs="i",
	        	main=sprintf("M=%1.2f", as.numeric(em)),
	        	colramp=colorRampPalette(c(addalpha('white', 0), addalpha(cols[i], 0.9)), alpha=T),
	        	nrpoints=0,
			add=T
		)
		points(dRats[i], bRats[i], col=cols[i], pch=19)
		#
		i = i+1
	}
	dev.off()
	#
	system('convert -delay 40 ./pics/M_0.04_0.30CloudGif*.png ./pics/M_0.04_0.30Cloud.gif')
}


##
#n = 10000
##
#l0 = cbind(rnorm(n, 10, 100), rnorm(n, 10, 100))
#l1 = cbind(rnorm(n, 20, 10), rnorm(10000, 20, 10))
#l2 = cbind(rnorm(n, -20, 10), rnorm(10000, 20, 10))
#l3 = cbind(rnorm(n, 20, 10), rnorm(10000, -20, 10))
#l4 = cbind(rnorm(n, -20, 10), rnorm(10000, -20, 10))

#
#PLOT
#

##
#dev.new()
#cols = brewer.pal(10, 'Spectral')
#smoothScatter(l0, 
#	xlim = c(min(l0[,1]), max(l0[,1])), xaxs="i",
#	ylim = c(min(l0[,2]), max(l0[,2])), yaxs="i",
#	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[1], 0.9)), alpha=T),
#	nrpoints = 0
#)
#smoothScatter(l1, 
#	add=T,
#	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[2], 0.9)), alpha=T),
#	nrpoints = 0
#)
#smoothScatter(l2, 
#	add=T,
#	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[3], 0.9)), alpha=T),
#	nrpoints = 0
#)
#smoothScatter(l3, 
#	add=T,
#	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[4], 0.9)), alpha=T),
#	nrpoints = 0
#)
#smoothScatter(l4, 
#	add=T,
#	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[5], 0.9)), alpha=T),
#	nrpoints = 0
#)

