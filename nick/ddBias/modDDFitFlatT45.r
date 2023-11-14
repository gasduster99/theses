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
FtFmsy = rep(1, TT)
##
#mid = round(TT/2)
#cMax = 2 #4*M
#cMin = 0.2 #M/4
#bb = log(cMin/cMax)/(tt-mid)
#aa = exp(log(cMax)-bb*mid)
#rSlope = (cMax-1)/(mid-1)
#rb = 2*rSlope*mid + cMax-rSlope*mid
#FtFmsy = (aa*exp(bb*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
##
#FtFmsy = c(FtFmsy, rep(1, 15))
#TT = length(FtFmsy)
#time = tt:TT

#DD MODEL STUFF

#A-0.5AS15K0.1
aS = 2	#2 #10 #0.1
a0 = -1		#-0.25 #-0.5 #-1   #-2
M  = 0.2
kappa = 0.10	#10 #0.1
WW = 1
ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
#
B0 = 10000
#
aMin = M*(M+kappa)/kappa/WW/(1+M*ww/kappa/WW)

#
#SIM
#

#a place to store data
##place = "./modsDDExpT45N300AS0.1K10N56/"	#zooid-3  #"./modsDDExpT45N150A-0.5AS15K0.1/"
#place = "./modsDDFlatT45N300A0-1AS10K0.1/"	#zooid-4  #"./modsDDExpT45N300AS10K0.1/" #"./modsDDExpT45N150A-1AS15K0.1/"
#place = "./modsDDFlatT45N300A0-1AS0.1K10/"  	#zooid-2  #"./modsDDExpT45N150A-1AS2/"     
##place = "./modsDDExpT45N300AS1K1N28/"    	#zooid-1  #"./modsDDExpT45N150A-0.5AS2/"   

##maybe I'll fit this later
#place = "./modsDDFlatT45N150A0-1AS0.1K10/"; rv=F;   #zooid1
place = "./modsDDFlatT45N150A0-1AS2K0.1/"; rv=T;   #zooid2
###most interesting start here
##place = "./modsDDFlatT45N150A0-1AS3K0.1/"; rv=F;   #zooid3
##place = "./modsDDFlatT45N300A0-1AS10K0.1N84/"; rv=T;   #zooid4

odeMethod = "lsode" #"radau" #

#
datFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("datGen*.rda")))

#
jobis = 1:length(datFiles)
#
job12 = seq(1, floor(length(datFiles)/2))
job22 = seq(length(datFiles), floor(length(datFiles)/2)+1)
jobsh = c(rbind(job12, job22))

#
registerDoParallel(46)
opts = list(preschedule=F)
foreach(i=rev(jobsh)*rv + jobsh*(1-rv), .options.multicore = opts) %dopar% {
#foreach(i=sort(1:length(datFiles), decreasing=rv), .options.multicore = opts) %dopar% {
#for(i in 1:length(datFiles)){
	#
        #DATA
        #

        #read in data from design locations
        datGen = readRDS(datFiles[i])
        datGen$M = M
	datGen$a0 = a0
	datGen$aS = aS
	datGen$kappa = kappa
        datGen$time = 1:TT
        datGen$catch = FtFmsy
        datGen$iterate(odeMethod)
        datGen$save( datFiles[i] )
	
	#
        #SKIP
        #

        #
        fileFit = gsub("datGen", "fit", datFiles[i])
        if( file.exists(fileFit) ){
                #writeLines(sprintf('\nSKIP Xi: %s, Zeta: %s\n', xiSims[j], zetaSims[i]))
                next
        }
	
	#
        #DATA
        #

        #make data to fit
        cpue = rlnorm(TT, datGen$lq+log(datGen$B), exp(datGen$lsdo))

	#
	#FIT
	#
	
	#
	lbetaBH = log(getBeta(B0, M, kappa, ww, WW, exp(datGen$lalpha), -1))
	fit = ddModel$new( derivs=der,
		N0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){#virgin numbers
		        #
		        alpha = exp(lalpha)
		        beta  = exp(lbeta)
		        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
		        #
		        BZero = BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
		        (alpha*BZero*( 1-beta*gamma*BZero )^(1/gamma))/M
		
		},
		B0Funk=function(lalpha, lbeta, gamma, M, WW, kappa, a0){#virgin biomass
		        #
		        alpha = exp(lalpha)
		        beta = exp(lbeta)
		        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
		        #
		        BBar(0, M, kappa, ww, WW, alpha, beta, gamma)
		},
		time=1:TT, catch=FtFmsy, aS=aS, a0=a0, M=M, WW=WW, kappa=kappa,#constants
		lalpha=log(exp(datGen$lalpha)), lbeta=lbetaBH, gamma=-1,         	#parameters
		lq=log(0.00049), lsdo=log(0.01160256),                  #nuisance parameters
		xi=datGen$xi, zeta=datGen$zeta                  	#other incidentals to carry along
	)
	fit$iterate(odeMethod)

	##
	#plot(cpue)
	#fit$plotMean(add=T)
	#fit$printSelf()

	
	#optimization   
        optAns = fit$optimize(cpue,
                c('lsdo', 'lalpha'),
                lower   = c(log(0.001), log(aMin)),
                upper   = c(log(1), log(10)),
                gaBoost = list(run=10, parallel=F, popSize=10^3), #10^4),
                persistFor = 5,
                fitQ    = F
        )
        
	##
	#fit$plotMean(add=T, col='red')
	#fit$printSelf()
			
	#
        optAns = fit$optimize(cpue,
                c('lsdo', 'lalpha', 'lbeta'),
                lower   = c(log(0.001), log(aMin), -10),
                upper   = c(log(1), log(10), -2),      #log(getBeta(100, -1, M, P0))
                gaBoost = list(run=100, parallel=F, popSize=10^3),#10^4),
                persistFor = 5,
                fitQ    = F
        )
        #get hessian if possible
        tryCatch({
                optAns = fit$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta'),
                        lower   = c(log(0.001), log(aMin), -10),
                        upper   = c(log(1), log(10), -2),
                        cov     = T,
                        fitQ    = F
                )
        }, error=function(err){
                writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
                optAns = fit$optimize(cpue,
                        c('lsdo', 'lalpha', 'lbeta'),
                        lower   = c(log(0.001), log(aMin), -10),
                        upper   = c(log(1), log(10), -2),
                        cov     = F,
                        fitQ    = F
                )
        })
	
	#convenience
        fit$alpha = exp(fit$lalpha)
        fit$beta  = exp(fit$lbeta)
        fit$sdo   = exp(fit$lsdo)
        fit$q     = exp(fit$lq)
        #save
        fit$save( fileFit )
        ##plot
        #fit$plotMean(add=T, col='blue')
        #fit$plotBand(col='blue', alpha=50)
}



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

