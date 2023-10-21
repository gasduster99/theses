rm(list=ls())

#
library(Ryacas)
library(pracma)
library(foreach)
library(parallel)
library(doParallel)

#
source('ddClass0.1.1.r')

#
#FUNCTIONS
#

#Equilibrium Equations

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
Bbar_r = as_r( with_value(Bbar, "F", ysym("FF")) )
BBar = function(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#
FBbar = ysym("F*B")
FBbar = with_value(FBbar, "B", Bbar)
FBbar = with_value(FBbar, "F", ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
dFBdF_r = as_r(dFBdF)
FDebug = function(FF, M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
	eval(dFBdFexpr)
}
FMsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10))$root
}

#beta does not matter for either of getAlphaFmsy or getGammaFmsy
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma, dFBdFexpr=dFBdF_r){
        #
        #capture.output(
        out <- tryCatch({
                uniroot(function(alpha){ eval(dFBdFexpr) }, c(eps(), 100))$root
        }, error=function(err){
                out = NA
        })#, file="/dev/null")
        #
        return(out)
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta, dFBdFexpr=dFBdF_r){
        #FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(dFBdFexpr) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
        f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
        uniroot(f, c(0, 10), tol=eps())$root
}

#gamma|alpha, zeta
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
        f = function(g){
                BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - zeta
        }
        uniroot(f, c(-10, 10))$root
}

#
getZeta = function(FF, M, k, w, W, alpha, beta, gamma){
        BBar(FF, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
}

#
vbGrow = function(a, k, W, a0){
	W*(1-exp(-k*(a-a0)))
}

#
getZetaBH = function(x, M, k, W, aS, a0){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
        #
        gamma = -1
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaBH = Vectorize(getZetaBH, "x")

#Differential Equation 

#
der = function(t, Y, lalpha, lbeta, gamma, aS, a0, WW, kappa, catch, B0){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        N = Y[1]
        B = Y[2]
        #
        if( (t-aS)<1){
                Blag = B0
        }else{
                Blag = lagvalue(t-aS)[2]
        }
        #
        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        alpha = exp(lalpha)
        beta = exp(lbeta)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
        #
        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
        FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
        #
        out = c(N=NA, B=NA)
        out[1] = R - (M+FF)*N
        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B
        #
        return( list(out) )
}

#
#INIT
#

#
TT = 31 #length(hake)
tt = 1
time = tt:TT
#
mid = round(TT/2)
cMax = 2 #4*M
cMin = 0.2 #M/4
bb = log(cMin/cMax)/(tt-mid)
aa = exp(log(cMax)-bb*mid)
rSlope = (cMax-1)/(mid-1)
rb = 2*rSlope*mid + cMax-rSlope*mid
FtFmsy = (aa*exp(bb*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
#
FtFmsy = c(FtFmsy, rep(1, 15))
TT = length(FtFmsy)
time = tt:TT

#DD MODEL STUFF

#
aS = 2
a0 = -0.5 #-0.25 #-0.5 #-1   #-2
M  = 0.2
kappa = 1 #0.2
WW = 1
ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
#
B0 = 10000

#
#SIM
#

#NOTE: should the server list fits? or datGen?
served = "./expFitters.csv"
#a place to store data
#place = "./modsDDExpT45N150A-0.5AS2/" #'./modsDDExpT45N150Wide/' #'./test/'#
odeMethod = "lsode" #"radau" #

#NOTE: should the server list fits? or datGen?
datFiles = unlist(read.csv(served, header=F))
#sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("datGen*.rda")))

##
#i = 126 #(3.19, 0.487) #66 yes (1.767, 0.48) #65 yes #15 no (0.577, 0.619) #67 no (1.79, 0.608)
#

#
jobis = 1:length(datFiles)
#
job12 = seq(1, floor(length(datFiles)/2))
job22 = seq(length(datFiles), floor(length(datFiles)/2)+1)
jobsh = c(rbind(job12, job22))

##
#registerDoParallel(8) #46)
#opts = list(preschedule=F)
#foreach(i=(1:length(datFiles)), .options.multicore = opts) %dopar% {
#foreach(i=rev(1:length(datFiles)), .options.multicore = opts) %dopar% {
for(i in (1:length(datFiles))[1]){
	#
        #DATA
        #

        #read in data from design locations
        datGen = readRDS(datFiles[i])
        #datGen$M = M
	#datGen$a0 = a0
	#datGen$aS = aS
	#datGen$kappa = kappa
        #datGen$time = 1:TT
        #datGen$catch = FtFmsy
        #datGen$iterate(odeMethod)
        #datGen$save( datFiles[i] )
	
	#
        #SKIP
        #

	#
	#read in fitBH and add fitBHKA and add fit3, fit3KA
	
	#if fit not found estimate fit, else read-in
		#fit fitBHKA from fit2
	#if fit3 not found estimate fit3, else read in 
		#fit fit3KA from fit3

	##
	#fileFit3 = gsub("datGen", "fit3", datFiles[i])
	#if( file.exists(fileFit3) ){
        #        #writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
        #        next
        #}
	
	##
        #fileFitBHKA = gsub("datGen", "fitBHKA", datFiles[i])
        #if( file.exists(fileFitBHKA) ){
        #        #writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
        #        next
        #}

	#
	fileFitBHKA = gsub("datGen", "fitBHKA", datFiles[i])
	fileFit3 = gsub("datGen", "fit3", datFiles[i])
	
	#
        fileFitBH = gsub("datGen", "fit", datFiles[i])
        if( !file.exists(fileFitBH) ){
                #writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
                next
        }
	
	#
	fileFit3KA = gsub("datGen", "fit3KA", datFiles[i])
	if( file.exists(fileFit3KA) ){
                #writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
                next
        }
	
	#
        #DATA
        #

        #make data to fit
        cpue = rlnorm(TT, datGen$lq+log(datGen$B), exp(datGen$lsdo)*10)
	#
        fitBH = readRDS(fileFitBH)


	#
	#FIT BHKA
	#
	
	#
	#fitBH = readRDS(fileFitBH)
	#lbetaBH = log(getBeta(B0, M, kappa, ww, WW, exp(datGen$lalpha), -1))
	fitBHKA = ddModel$new( derivs=der,
		N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin numbers
		        #
		        alpha = exp(lalpha)
		        beta  = exp(lbeta)
		        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
		        #
		        BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
		        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
		
		},
		B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin biomass
		        #
		        alpha = exp(lalpha)
		        beta = exp(lbeta)
		        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
		        #
		        BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
		},
		time=fitBH$time, catch=fitBH$catch, a0=fitBH$a0, M=fitBH$M, 	#constants 
		aS=fitBH$aS, WW=fitBH$WW, kappa=fitBH$kappa,			#growth
		lalpha=fitBH$lalpha, lbeta=fitBH$lbeta, gamma=-1,         	#recruitment
		lq=fitBH$lq, lsdo=fitBH$lsdo,                  			#nuisance
		xi=fitBH$xi, zeta=fitBH$zeta, cpue=cpue              	#other incidentals to carry along
	)
	fitBHKA$iterate(odeMethod)	
	
	#
        optAns = fitBHKA$optimize(cpue,
                c('lsdo', 'lalpha', 'lbeta', 'kappa', 'aS'),
                lower   = c(log(0.001), log(M+0.1), -10, 0.1, 0.1),
                upper   = c(log(1), log(10), -2, 100, 100),       #log(getBeta(100, -1, M, P0))
                #gaBoost = list(run=10, parallel=T, popSize=10^3), 
                #gaBoost = list(run=10, parallel=F, popSize=10^3),#10^4),
                #persistFor = 5,
                fitQ    = F
        )
	#get hessian if possible
        tryCatch({
                optAns = fitBHKA$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta', 'kappa', 'aS'),
                        lower   = c(log(0.001), log(M+0.1), -10, 0.1, 0.1),
                        upper   = c(log(1), log(10), -2, 100, 100),
                        cov     = T,
                        fitQ    = F
                )
        }, error=function(err){
                writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
                optAns = fitBHKA$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta', 'kappa', 'aS'),
                        lower   = c(log(0.001), log(M+0.1), -10, 0.1, 0.1),
                        upper   = c(log(1), log(10), -2, 100, 100),
                        cov     = F,
                        fitQ    = F
                )
        })
	
	##
	#fit$plotMean(add=T, col='blue')
	#fit$printSelf()

	#convenience
	fitBHKA$alpha = exp(fitBHKA$lalpha)
        fitBHKA$beta  = exp(fitBHKA$lbeta)
        fitBHKA$sdo   = exp(fitBHKA$lsdo)
        fitBHKA$q     = exp(fitBHKA$lq)
        #save
        fitBHKA$save( fileFitBHKA )
        
	#
	#FIT3 
	#

	#
	wwDat = vbGrow(datGen$aS, datGen$kappa, datGen$WW, datGen$a0)
	lbeta3 = log(getBeta(datGen$B0, datGen$M, datGen$kappa, wwDat, datGen$WW, exp(datGen$lalpha), datGen$gamma))
	fit3 = ddModel$new( derivs=der,
		N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin numbers
		        #
		        alpha = exp(lalpha)
		        beta  = exp(lbeta)
		        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
		        #
		        BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
		        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
		
		},
		B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin biomass
		        #
		        alpha = exp(lalpha)
		        beta = exp(lbeta)
		        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
		        #
		        BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
		},
		time=fitBH$time, catch=fitBH$catch, a0=fitBH$a0, M=fitBH$M, 	#constants 
		aS=fitBH$aS, WW=fitBH$WW, kappa=fitBH$kappa,			#growth
		lalpha=fitBH$lalpha, lbeta=lbeta3, gamma=datGen$gamma,  #-1, 	#recruitment
		lq=fitBH$lq, lsdo=fitBH$lsdo,                  			#nuisance
		xi=fitBH$xi, zeta=fitBH$zeta, cpue=cpue               	#other incidentals to carry along
	)
	fit3$iterate(odeMethod)

	###
        ##optAns = fit3$optimize(cpue,
        ##        c('gamma'),
        ##        lower   = c(-2),
        ##        upper   = c(2),      #log(getBeta(100, -1, M, P0))
        ##        gaBoost = list(run=10, parallel=F, popSize=10^3),#10^4),
        ##        persistFor = 5,
        ##        fitQ    = F
        ##)

	##
        #optAns = fit3$optimize(cpue,
        #        c('lsdo', 'lalpha'),
        #        lower   = c(log(0.001), log(M+0.1)),
        #        upper   = c(log(1), log(10)),      #log(getBeta(100, -1, M, P0))
        #        gaBoost = list(run=10, parallel=T, popSize=10^3),#10^4),
	#	#gaBoost = list(run=10, parallel=F, popSize=10^3),#10^4),
        #        persistFor = 5,
        #        fitQ    = F
        #)

	##
        #optAns = fit3$optimize(cpue,
        #        c('lsdo', 'lalpha', 'lbeta', 'gamma'),
        #        lower   = c(log(0.001), log(M+0.1), -10, -2),
        #        upper   = c(log(1), log(10), -2, 2),      #log(getBeta(100, -1, M, P0))
        #        gaBoost = list(run=10, parallel=T, popSize=10^3),#10^4),
	#	#gaBoost = list(run=10, parallel=F, popSize=10^3),#10^4),
        #        persistFor = 5,
        #        fitQ    = F
        #)
	##get hessian if possible
        #tryCatch({
        #        optAns = fit3$optimize(cpue,
        #                c('lsdo', 'lalpha', 'lbeta', 'gamma'),
        #                lower   = c(log(0.001), log(M+0.1), -10, -2),
        #                upper   = c(log(1), log(10), -2, 2),
        #                cov     = T,
        #                fitQ    = F
        #        )
        #}, error=function(err){
        #        writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
        #        optAns = fit3$optimize(cpue,
        #                c('lsdo', 'lalpha', 'lbeta', 'gamma'),
        #                lower   = c(log(0.001), log(M+0.1), -10, -2),
        #                upper   = c(log(1), log(10), -2, 2),
        #                cov     = F,
        #                fitQ    = F
        #        )
        #})
	#
	###
        ##fit$plotMean(add=T, col='blue')
        ##fit$printSelf()

        ##convenience
        #fit3$alpha = exp(fit3$lalpha)
        #fit3$beta  = exp(fit3$lbeta)
        #fit3$sdo   = exp(fit3$lsdo)
        #fit3$q     = exp(fit3$lq)
        ##save
        #fit3$save( fileFit3 )
        ###plot
        ##fit$plotMean(add=T, col='blue')
        ##fit$plotBand(col='blue', alpha=50)

	##
	##FIT 3KA
	##
	#
	##
	##fitBH = readRDS(fileFitBH)
	##lbetaBH = log(getBeta(B0, M, kappa, ww, WW, exp(datGen$lalpha), -1))
	#fit3KA = ddModel$new( derivs=der,
	#	N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin numbers
	#	        #
	#	        alpha = exp(lalpha)
	#	        beta  = exp(lbeta)
	#	        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
	#	        #
	#	        BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
	#	        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
	#	
	#	},
	#	B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0, aS){#virgin biomass
	#	        #
	#	        alpha = exp(lalpha)
	#	        beta = exp(lbeta)
	#	        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
	#	        #
	#	        BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
	#	},
	#	time=fit3$time, catch=fit3$catch, a0=fit3$a0, M=fit3$M, 	#constants 
	#	aS=fit3$aS, WW=fit3$WW, kappa=fit3$kappa,			#growth
	#	lalpha=fit3$lalpha, lbeta=fit3$lbeta, gamma=fit3$gamma,         	#recruitment
	#	lq=fit3$lq, lsdo=fit3$lsdo,                  			#nuisance
	#	xi=fit3$xi, zeta=fit3$zeta, cpue=cpue               	#other incidentals to carry along
	#)
	#fit3KA$iterate(odeMethod)	
	fit3KA = fit3
	
	#
        optAns = fit3KA$optimize(cpue,
                c('kappa', 'aS'),
                lower   = c(0.1, 0.1),
                upper   = c(100, 100),      #log(getBeta(100, -1, M, P0))
                gaBoost = list(run=10, parallel=T, popSize=10^3),#10^4),
		#gaBoost = list(run=10, parallel=F, popSize=10^3),#10^4),
                persistFor = 5,
                fitQ    = F
        )
	
	#
	plot(cpue)
	fit3KA$plotMean(add=T, col='cyan')
	fit3KA$printSelf()
	
	#
        optAns = fit3KA$optimize(cpue,
                c('lsdo', 'lalpha', 'lbeta', 'gamma', 'kappa'),
                lower   = c(log(0.001), log(M+0.1), -10, -2, 0.1),
                upper   = c(log(1), log(10), -2, 2, 100),      #log(getBeta(100, -1, M, P0))
                gaBoost = list(run=10, parallel=T, popSize=10^3),#10^4),
		#gaBoost = list(run=10, parallel=F, popSize=10^3),#10^4),
                persistFor = 5,
                fitQ    = F
        )
	#get hessian if possible
        tryCatch({
                optAns = fit3KA$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta', 'gamma', 'kappa', 'aS'),
                        lower   = c(log(0.001), log(M+0.1), -10, -2, 0.1, 0.1),
                        upper   = c(log(1), log(10), -2, 2, 100, 100),
                        cov     = T,
                        fitQ    = F
                )
        }, error=function(err){
                writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
                optAns = fit3KA$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta', 'gamma', 'kappa', 'aS'),
                        lower   = c(log(0.001), log(M+0.1), -10, -2, 0.1, 0.1),
                        upper   = c(log(1), log(10), -2, 2, 100, 100),
                        cov     = F,
                        fitQ    = F
                )
        })
	
	###
	##fit$plotMean(add=T, col='blue')
	##fit$printSelf()

	##convenience
	#fit3KA$alpha = exp(fit3KA$lalpha)
        #fit3KA$beta  = exp(fit3KA$lbeta)
        #fit3KA$sdo   = exp(fit3KA$lsdo)
        #fit3KA$q     = exp(fit3KA$lq)
        ##save
        #fit3KA$save( fileFit3KA ) 
}




##
	#plot(cpue)
	#fit$plotMean(add=T)
	#fit$printSelf()

	#
	##optimization   
        #optAns = fit$optimize(cpue,
        #        c('lsdo', 'lalpha'),
        #        lower   = c(log(0.001), log(M+kappa)),
        #        upper   = c(log(1), log(100)),
        #        gaBoost = list(run=10, parallel=T, popSize=10^3), #10^4),
        #        persistFor = 5,
        #        fitQ    = F
        #)
        #
	##
	#fit$plotMean(add=T, col='red')
	#fit$printSelf()
	#		
	##
        #optAns = fit$optimize(cpue,
        #        c('lsdo', 'lalpha', 'lbeta'),
        #        lower   = c(log(0.001), log(M+kappa), -10),
        #        upper   = c(log(1), log(100), -2),      #log(getBeta(100, -1, M, P0))
        #        gaBoost = list(run=100, parallel=T, popSize=10^3),#10^4),
        #        persistFor = 5,
        #        fitQ    = F
        #)
        ##get hessian if possible
        #tryCatch({
        #        optAns = fit$optimize(cpue,
        #                c('lsdo', 'lalpha', 'lbeta'),
        #                lower   = c(log(0.001), log(M+kappa), -10),
        #                upper   = c(log(1), log(100), -2),
        #                cov     = T,
        #                fitQ    = F
        #        )
        #}, error=function(err){
        #        writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
        #        optAns = fit$optimize(cpue,
        #                c('lsdo', 'lalpha', 'lbeta'),
        #                lower   = c(log(0.001), log(M+kappa), -10),
        #                upper   = c(log(1), log(100), -2),
        #                cov     = F,
        #                fitQ    = F
        #        )
        #})

	##
	#fit$plotMean(add=T, col='green')
	#fit$printSelf()


#fit = prodModel$new(
#                dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)}, #model
#                time=1:TT, catch=FtFmsy, M=M,                                #BH        #constants
#                alpha=datGen$alpha, beta=getBeta(datGen$alpha, -1, M, P0), gamma=-1,    #parameters
#                lalpha=datGen$lalpha, lbeta=log(getBeta(datGen$alpha, -1, M, P0)),      #reparameterize
#                lq=log(0.00049), lsdo=log(0.01160256),          #nuisance parameters
#                xi=datGen$xi, zeta=datGen$zeta                  #other incidentals to carry along
#        )



#
#	#
#	lalphaBH2 = log(getAlphaFmsy(datGen$xi*M, M, kappa, ww, WW, 1, -1))
#	lbetaBH2 = log(getBeta(B0, M, kappa, ww, WW, exp(lalphaBH2), -1))
#	fit2 = ddModel$new( derivs=der,
#		N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){#virgin numbers
#		        #
#		        alpha = exp(lalpha)
#		        beta  = exp(lbeta)
#		        ww = WW*(1-exp(-kappa*a0))
#		        #
#		        BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
#		        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
#		
#		},
#		B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){#virgin biomass
#		        #
#		        alpha = exp(lalpha)
#		        beta = exp(lbeta)
#		        ww = WW*(1-exp(-kappa*a0))
#		        #
#		        BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
#		},
#		time=1:TT, catch=FtFmsy, a0=a0, M=M, WW=WW, kappa=kappa,#constants
#		lalpha=lalphaBH2, lbeta=lbetaBH2, gamma=-1,         	#parameters
#		lq=log(0.00049), lsdo=log(0.01160256),                  #nuisance parameters
#		xi=datGen$xi, zeta=datGen$zeta                  	#other incidentals to carry along
#	)
#	fit2$iterate(odeMethod)
#
#
#
#
#

