rm(list=ls())

#
library(tgp)
#library(parallel)
library(bigmemory)
library(microbenchmark)
suppressMessages(library(foreach, quietly=FALSE))
suppressMessages(library(doParallel, quietly=FALSE))

#
#FUNCTIONS
#

#
michal = function(xx, m=10, p=2){
        ##########################################################################
        #
        # INPUTS:
        #
        # xx = c(x1, x2)
        # m = constant (optional), with default value 10
        #
        ##########################################################################

        #
        xx = matrix(xx, ncol=p)
        ii = 1:p
        sum = rowSums(sin(xx) * (sin(sweep(xx^2/pi, 2, ii, '*')))^(2*m))
        #
        y = -sum
        return(y)
}

#
grlee12 = function(xx){
        #xMin: 0.5485634
	#zMin: -0.8690111
	
	#
        xx = matrix(xx, ncol=1)
        #
        term1 = sin(10*pi*xx) / (2*xx)
        term2 = (xx-1)^4
        #
        y = term1 + term2
        return(y)
}

#
rosenbrock = function(x){
        #
        x = matrix(x, ncol=2)
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

#
rastrigin = function(x, p=2){
        #
        x = matrix(x, ncol=p)
        out = matrix(NaN, nrow=dim(x)[1], ncol=1)
        for (i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = 10*length(ex)+sum(ex^2-10*cos(2*pi*ex))
        }
        return(out)
}

#
levy = function(xx, d=2){
        ##########################################################################
        #
        # INPUT:
        #
        # xx = matrix(c(x1, x2, ..., xd), ncol=d)
        #
        ##########################################################################

        #
        xx = matrix(xx, ncol=d)
        #
        w = matrix(apply( xx, 2, function(x){1+(x-1)/4} ), ncol=d)
        term1 = (sin(pi*w[,1]))^2
        term3 = (w[,d]-1)^2 * (1+1*(sin(2*pi*w[,d]))^2)
        #
        wi  = matrix(w[,1:(d-1)], ncol=d-1)
        sum = rowSums((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
        # 
        y = term1 + sum + term3
        #
        return(y)
}

#
optimStep = function(f, rect, trace=F){
	#
	xInitPerVol = 2
        xInit = xInitPerVol*prod(rect[,2]-rect[,1])	
	X = lhs(xInit, rect)
	tgp::optim.step.tgp( f, X=X, Z=f(X), rect=rect, improv=c(1,1), verb=0, NN=200, trace=trace )
}

#
#BENCH
#

#
fOnly = microbenchmark(
	"grlee12"    = grlee12(0.5485634),
	"rosenbrock" = rosenbrock(c(1, 1)),
	"michal"     = michal(c(2.20, 1.57)),
	"rastrigin"  = rastrigin(c(0, 0)),
	"levy"	     = levy(c(1, 1)),
	times        = 10^4
)

##
#fTGP = microbenchmark(
#        "grlee12"    = optimStep(grlee12, cbind(c(0.5), c(2.5))),
#        "rosenbrock" = optimStep(rosenbrock, cbind(c(-2, -3), c(2, 5))),
#        "michal"     = optimStep(michal, cbind(c(0, 0), c(pi, pi))),
#        "rastrigin"  = optimStep(rastrigin, cbind(c(-2.5, -2.5), c(2.5, 2.5))),
#        "levy"       = optimStep(levy, cbind(c(-10, -10), c(10, 10))),
#        times        = 10 
#)

#
names = c("grlee12", "rosenbrock", "michal", "rastrigin", "levy")
fs = list(
	grlee12=grlee12,
	rosenbrock=rosenbrock, 
	michal=michal, 
	rastrigin=rastrigin,
	levy=levy
)
rects = list(
	grlee12=cbind(c(0.5), c(2.5)),
	rosenbrock=cbind(c(-2, -3), c(2, 5)),
	michal=cbind(c(0, 0), c(pi, pi)),
	rastrigin=cbind(c(-2.5, -2.5), c(2.5, 2.5)),
	levy=cbind(c(-10, -10), c(10, 10))
)
vols = sapply(rects, function(r){prod(r[,1]-r[,2])})

#
trace = F
#
threads = min(8, length(names))
#
pListDef = big.matrix(threads, 1,
        init           = -1,
        backingfile    = 'pList.bin',
        descriptorfile = "pList.desc"
)

#
M = 8
threads = 8
#
registerDoParallel(cores=threads)
out = foreach( m=1:M )%dopar%{
	#parallel
        pList = attach.big.matrix("pList.desc")
        pid = Sys.getpid()
        who = which(as.matrix(pList)==-1)
        if( m<=threads ){ pList[m] = pid
        } else{ pList[who[sample(1:length(who),1)]] = pid }
        rank = which(as.matrix(pList)==pid)-1
	#move into a clean work space for TGP
	outPath = getwd()
        dir.create(sprintf('%s/rank%02d', outPath, rank))
        setwd(sprintf('%s/rank%02d', outPath, rank))
	#benchmark
	#xInitPerVol = 2
	#xInit = xInitPerVol*prod(rects[[n]][,2]-rects[[n]][,1])
	fOut = microbenchmark(
		optimStep(fs[[1]], rects[[1]], trace), 
		optimStep(fs[[2]], rects[[2]], trace),
		optimStep(fs[[3]], rects[[3]], trace),
		optimStep(fs[[4]], rects[[4]], trace),
		#optimStep(fs[[5]], rects[[5]], trace),
		unit='s', 
		times=1
	)
	#
	return( list(fOut=summary(fOut)) )
}
##
#fOut = sapply(out, function(p){p['fOut']})
#fOut = do.call(rbind, fOut)
#fOut = fOut[,-1]
#rownames(fOut) = names
##
#xInit = sapply(out, function(p){p['xInit']})
#xInit = do.call(rbind, xInit)
##
#fOut = cbind(fOut, xInit)







##
#ffOut = c()
#for(n in names){
#	out = microbenchmark(optimStep(fs[[n]], rects[[n]]), times=2)
#	ffOut = rbind(ffOut, summary(out))
#}
