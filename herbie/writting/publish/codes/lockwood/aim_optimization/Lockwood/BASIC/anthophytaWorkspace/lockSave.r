rm(list=ls())

library(tgp)

source('tSource.r')
source('alerts.r')

#
runlockBasic2 = function(rates=c(rep(10000,2))){ 
	#
	vLoss = 40000
	#vector input version
  	write(c(6,rates,10000,10000,10000,10000),file="input.tst",ncol=1)
  	system("./RunLock input.tst output.tst 1")
  	output = scan(file="output.tst")
  	Y = output[1] + vLoss*(abs(output[2])>0) + vLoss*(abs(output[3])>0) + 2*abs(output[2]) + 2*abs(output[3])
  	#
	return( Y )
}

#
runlockBasic6 = function(rates=c(rep(10000,6))){ 
	#
	vLoss = 40000
	#vector input version
	write(c(6, rates),file="input.tst",ncol=1)
  	system("./RunLock input.tst output.tst 1")
  	output = scan(file="output.tst")
  	Y = output[1] + vLoss*(abs(output[2])>0) + vLoss*(abs(output[3])>0) + 2*abs(output[2]) + 2*abs(output[3])
  	#
	return(Y)

}

#runlockBasic6 = function(rates=c(rep(10000,6))){
#	write(c(6,rates),file="input.tst",ncol=1)
# 	system("./RunLock input.tst output.tst 1")
#  	output = scan(file="output.tst")
#	#
#  	return( data.frame(cost=output[1],plume.a=output[2],plume.b=output[3]) )
#}

#typically keep will W*sLen (ie. number of points looking back (W) by the samples per point)
slideWindow = function(oldStuff, newStuff, keep){
        out = c(oldStuff, newStuff)
        outl = length(out)
        left = max((outl-(keep-1)), 1)
        out = out[left:outl]

        return( out )
}

#
#

#
f = runlockBasic2                #rosenbrock#rastrigin#easom#shekel#
rect = cbind(c(-2, -3), c(2, 5)) #cbind(c(-2, -3), c(2, 5))

#
em = dim(rect)[1]
#sample size per iteration
n = 1
#intially sampling the points 
X = lhs(40, rect)
#making the intial "sample" by evaluating the function at the sample points
Z = f(X)
#sweeps of optimization
M = 60#150#120
#maximum window size
W = 10
#signifigance level
sig = 3.

#intializing
Zmax = c(min(Z))
wSamples = c()
samples = c()
sLen = NULL
maxes = matrix(NaN, nrow=M, ncol=1)
out = NULL
#
dev.new(width=14, height=14/3)
#
meat = function(init, it){
        #[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:mcmcSize 
        X = init[[1]]
        Z = init[[2]]
        Zmax = init[[3]]
        samples = init[[4]]
        maxes = init[[5]]
        out = init[[6]]
        sLen = init[[7]]
        #               
        out = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=out, improv=c(1,n), trace=T, verb=0 )
        ex = matrix(out$X, ncol=em)
        fex = f(ex)
        X = rbind(X, ex)
        Z = c(Z, fex)
        Zmax = c(Zmax, min(Z))
        #
        improvSamples = out$obj$trace$preds$improv
        EimprovAll = out$obj$improv
        maxI = which( EimprovAll$rank==1 )
        maxSamples = unlist( improvSamples[maxI] )
        m = mean(maxSamples)
        maxes[it] = m
        #
        sLen = length(maxSamples)
        wSamples = slideWindow(wSamples, maxSamples, W*sLen)
        samples = c(samples, maxSamples)
        #
        wm = mean(wSamples)
        #
        samLower = quantile(wSamples, probs=pnorm(-sig, 0, 1))
        samUpper = quantile(wSamples, probs=pnorm(sig, 0, 1))
        #
        #OUTPUT
        layout( matrix(c(1, 2, 3), nrow=1, ncol=3) )
        #
        #Figure 1
        plot( seq(1, it+1), Zmax,
                "l",
                xlab="Sweeps",
                ylab="Value",
                main="Best Z Value"
        )
	#
        #Figure 2
        L = length(seq(1, it))
        left = max((L-(W-1)), 1)
        plot( seq(1, it), maxes[1:it],
                main=sprintf('Sweeps: %d\n', it),
                ylab='Max Improvement\n',
                col=c(rep("black", L-1), "green"),
                ylim=c(0, max(maxes[1:it], samUpper))
        )
        #Sample Intervals
        segments(    1, samUpper, left, samUpper, col="blue", lty=2)
        segments(    1, samLower, left, samLower, col="blue", lty=2)
        segments( left, samUpper,   it, samUpper, col="blue")
        segments( left, samLower,   it, samLower, col="blue")
        #BarBar
        segments( left, wm, it, wm, col="black")
        #
        #Figure 3
        hist( wSamples, add=F, freq=F,
                main="Window Samples"
        )
        abline( v=samUpper, col="blue")
        abline( v=samLower, col="blue")
        abline( v=m, col='green' )
        abline( v=wm, col='black' )
        legend("topright",
                legend=c(sprintf("%0.4f%s (%d\\sigma)\nLog-N Interval\n", pnorm(sig, 0, 1)-pnorm(-sig, 0, 1), "%", sig), sprintf("%0.4f%s (%d\\sigma)\nSample Interval\n", pnorm(sig, 0, 1)-pnorm(-sig, 0, 1), "%",sig), "x-bar", "x-bar-bar"), 
                col=c("red", "blue", "green", "black"),
                lty=1
        )
        #
        init = list(X, Z, Zmax, samples, maxes, out, sLen)
}

#
ZMaxSea = list(Zmax)
samplesSea = list(samples)
maxesSea = list(maxes)
sLenSea = list(sLen)
#[[1]][[i]]:ZMaxSea; [[2]][[i]]:samplesSea; [[3]][[i]]:maxesSea; [[4]][[i]]:sLen:
seaSave = list(ZMaxSea, samplesSea, maxesSea, sLenSea)

#[[1]]:X; [[2]]:Z; [[3]]:Zmax; [[4]]:wSamples; [[5]]:maxes; [[6]]:out; [[7]]:sLen
init = list(X, Z, Zmax, samples, maxes, out, sLen)
for(it in seq(1, M)){
        init = meat(init, it)
        #[[1]][[i]]:ZMaxSea; [[2]][[i]]:samplesSea; [[3]][[i]]:maxesSea; [[4]][[i]]:sLen;       
        seaSave[[1]] = c( seaSave[[1]], list(init[[3]]) )
        seaSave[[2]] = c( seaSave[[2]], list(init[[4]]) )
        seaSave[[3]] = c( seaSave[[3]], list(init[[5]]) )
        seaSave[[4]] = c( seaSave[[4]], list(init[[7]]) )
}

what = "Lock2"  #"Rast" #"Easom"  #"Shekel" #"Rast"#"
size = "small"  #"Hard" #"Med"    #"Med10"  #"Hard2"#"
save.image(sprintf("seaSave%s%s.RData", what, size))

ding(1)


