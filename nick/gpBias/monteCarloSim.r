rm(list=ls())

#
library(pracma)
library(foreach)
library(parallel)
library(doParallel)
#library(mvtnorm)
#library(rootSolve)
#library(latex2exp)
#library(matrixStats)
#library(RColorBrewer)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
getFits = function(dir, xi, zeta, binTrk=2){
        #
        globber = sprintf("fit_xi%s_zeta%s_*.rda", round(xi, binTrk), round(zeta, binTrk))
        #
        sprintf( "%s%s", dir, list.files(path=dir, pattern=glob2rx(globber)) )
}

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
        if( gamma==-1 ){ Fmsy = sqrt(exp(lalpha)*M)-M
        } else {
                #
                FUpper = abs(exp(lalpha)-M)
                Fmsy = uniroot.all(function(ff){ a(gamma, ff, M)-exp(lalpha) }, c(0, FUpper))[1]
        }
        #
        R = exp(lalpha)*P*(1-exp(lbeta)*gamma*P)^(1/gamma)
        out = R - (M+Fmsy*C)*P
        #
        return( list(out) )
}

#
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) }

#
getBeta = function(alpha, gamma, M, B0){
        (1-(M/alpha)^gamma)/B0/gamma
}
b = Vectorize(getBeta, c("alpha", "gamma"))

#
#FRONT MATTER
#

#
mod = "FlatT30N150Wide" #"FlatT30N150WideExpHotStart" #
dirIn = sprintf("./modsSchnute%s/", mod)
#
datFiles = sprintf("%s%s", dirIn, list.files(path=dirIn, pattern=glob2rx("datGen*.rda")))
#
xis = unlist(sapply(datFiles, function(fn){
                dat = readRDS(fn)
                if( length(dat$zeta)!=0 & length(dat$xi)!=0){
                        return(dat$xi)
                }
        })
)
zetas = unlist(sapply(datFiles, function(fn){
                dat = readRDS(fn)
                if( length(dat$zeta)!=0 & length(dat$xi)!=0){
                        return(dat$zeta)
                }
        })
)

#
xiPoint = 3.4	#0.5   #(3.5-0.5)/2 + 0.5
zetaPoint = 0.4	#0.5 #(0.7-0.2)/2 + 0.2

#minimize (max, max) norm
norms = sqrt((xiPoint-xis)^2 + (zetaPoint-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]

#safely create dirOut
start = 1
dirOut = sprintf("./monteCarloHHard%s/", mod)
if( dir.exists(dirOut) ){
        #
        isGo = readline(sprintf("Overwrite %s Simulation?\nCtrl-C to Escape, Enter to Overwrite, or A to Append: ", dirOut))
        #
        if( toupper(isGo)=="A" ){
                #
                dir.create(dirOut, showWarnings=FALSE)
		#
		fits = getFits(dirOut, xi, zeta)
		if(length(fits)==0){ start=1 
		}else{
			#
			fits = gsub('.rda', '', fits) 
			fits = strsplit(fits, 'sim')
			fits = as.numeric(unlist(fits)[c(F,T)])
			start = max(fits)
		}
        }else{
                #
                unlink(dirOut, recursive=TRUE)
                dir.create(dirOut, showWarnings=FALSE)
        }
} else{ dir.create(dirOut, showWarnings=FALSE) }
#a filename template for the fitting model
fitFile = names(who)
fitFile = sub("datGen", "fit", fitFile)
fitFile = sub(dirIn, dirOut, fitFile)

#
datGen    = readRDS(names(who)) 
odeMethod = datGen$ODE_method
#
TT = length(datGen$time)
M  = datGen$M
P0 = datGen$N0
#
FtFmsy = rep(1, TT)

#
threads = 46
MM = 100 #threads
registerDoParallel(threads)
opts = list(preschedule=F)
foreach(i=rev(start:(start+MM)), .options.multicore = opts) %dopar% {
	#
	cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
	
	#
	fit = prodModel$new(
	        dNdt=dPdt, N0Funk=function(lalpha, lbeta, gamma, M){PBar(exp(lalpha), exp(lbeta), gamma, 0, M)},
	        time=1:TT, catch=FtFmsy, M=M,                                #BH   
	        alpha=datGen$alpha, beta=getBeta(datGen$alpha, -1, M, P0), gamma=-1,
	        lalpha=datGen$lalpha, lbeta=log(getBeta(datGen$alpha, -1, M, P0)),
	        lq=log(0.00049), lsdo=log(0.01160256),          #nuisance parameter
	        xi=datGen$xi, zeta=datGen$zeta                  #other incidentals 
	)
	fit$iterate(odeMethod)
	
	#optimization   
	optAns = fit$optimize(cpue,
	        c('lsdo', 'lalpha'),
	        lower   = c(log(0.001), log(M)),
	        upper   = c(log(1), log(100)),
	        gaBoost = list(run=10, parallel=F, popSize=10^4),
	        persistFor = 5,
	        fitQ    = F
	)
	#
	optAns = fit$optimize(cpue,
	        c('lsdo', 'lalpha', 'lbeta'),
	        lower   = c(log(0.001), log(M), -10),
	        upper   = c(log(1), log(100), -2),      #log(getBeta(100, -1, M, P0))
	        gaBoost = list(run=100, parallel=F, popSize=10^4),
	        persistFor = 5,
	        fitQ    = F
	)
	#get hessian if possible
	tryCatch({
	        optAns = fit$optimize(cpue,
	                c('lsdo', 'lalpha', 'lbeta'),
	                lower   = c(log(0.001), log(M), -10),
	                upper   = c(log(1), log(100), -2),
	                cov     = T,
	                fitQ    = F
	        )
	}, error=function(err){
	        writeLines( sprintf("\nNO HESSIAN AT xi: %s | zeta:%s", datGen$xi, datGen$zeta) )
	        optAns = fit$optimize(cpue,
	                c('lsdo', 'lalpha', 'lbeta'),
	                lower   = c(log(0.001), log(M), -10),
	                upper   = c(log(1), log(100), -2),
	                cov     = F,
	                fitQ    = F
	        )
	})
	#convenience
	fit$alpha = exp(fit$lalpha)
	fit$beta  = exp(fit$lbeta)
	fit$sdo   = exp(fit$lsdo)
	fit$q     = exp(fit$lq)
	fit$cpue  = cpue
	#save
	sFile = sprintf('_sim%s.rda', i)
	sFile = sub(".rda", sFile, fitFile)
        fit$save( sFile )
}



#plot(cpue)
#for(i in 1:50){
#	cpue = rlnorm(TT, datGen$lq+log(datGen$N), exp(datGen$lsdo))
#	points(cpue)
#}



