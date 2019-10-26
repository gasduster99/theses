#rm( list=ls() )

library(tgp)

source('tSource.r')
source('alerts.r')

rosenbrock = function(x){
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

rastrigin <- function(x){
        out = matrix(NaN, nrow=dim(x)[1], ncol=1)
        for (i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = 10*length(ex)+sum(ex^2-10*cos(2*pi*ex))
        }   
        return(out)
}


#give xBarBar as the repeated sampling mean
#give xDist as the MCMC predictive distribution of x
#give sigma as the level of control you desire as an even number (default sigma=6)
#return true if you pass test
test1 = function(xDist, xBarBar, sigma=6){
	sigma = sigma/2
	xBar = mean(xDist)
	s2 = var(xDist)
	n = length(xDist)
	SE = sqrt(s2)/sqrt(n)
	
	#a row for every control level
	#column one: LCL, column two: UCL 
	CL = matrix(NaN, nrow=sigma, ncol=2)
	for( i in seq(1, sigma) ){
		CL[i,] = c(xBarBar-i*SE, xBarBar+i*SE)
	}

	print( xBar )	
	print( CL[3,] )
	
	#greater than the upper limit
	#or less than the lower limit
	#negate to pass(return T) test if you fail either above condition 
	return( !(xBar>CL[sigma,2] || xBar<CL[sigma,1]) )
}

#typically keep will W*sLen (ie. number of points looking back (W) by the samples per point)
slideWindow = function(oldStuff, newStuff, keep){
	out = c(oldStuff, newStuff)
	outl = length(out)
	left = max((outl-(keep-1)), 1)
	out = out[left:outl]
	
	return( out )	
}


f = rastrigin #rosenbrock #

#domain boundaries
rect = cbind(c(-2, -2), c(2, 2))#cbind(c(-1,-1), c(1, 1))
em = dim(rect)[1]

#initially I use 40
#sample size per iteration
n = 1
#intially sampling the points 
X = lhs(40, rect)
#making the intial "sample" by evaluating the function at the sample points
Z = f(X)
#sweeps of optimization
M = 60
W = 10

dev.new(width=14, height=14/3)

Zmax = c(min(Z))
wSamples = c()
maxes = matrix(NaN, nrow=M, ncol=1)
maxesLog = matrix(NaN, nrow=M, ncol=1)
out = NULL
for( i in seq(1, M) ){
	#readline("Press any key to go to next.\n")
	out = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=out, improv=c(1,n), trace=T, verb=0 )
        ex = matrix(out$X, ncol=em)
	fex = f(ex)
        X = rbind(X, ex)     
        Z = c(Z, fex)
	Zmax = c(Zmax, min(Z))
        #prog = out$progress
        #improv = c(improv, prog$improv)	

	improvSamples = out$obj$trace$preds$improv
	EimprovAll = out$obj$improv
	maxI = which( EimprovAll$rank==1 )

	maxSamples = unlist( improvSamples[maxI] )
	m = mean(maxSamples)
	v = var(maxSamples)	
	mu = log( (m^2)/(v+m^2)^0.5 )
	s2 = log( 1 + (v/(m^2)) )
	upper = exp( mu + (1.96*(s2^0.5)) )
	lower = exp( mu - (1.96*(s2^0.5)) )
	
	#maxSamplesLog = log(maxSamples)
	maxes[i] = m #mean(maxSamples)#EimprovAll[maxI,1]
	#maxesLog[i] = mean(maxSamplesLog)

        sLen = length(maxSamples)
        wSamples = slideWindow(wSamples, maxSamples, W*sLen)
	
	wm = mean(wSamples)
	wv = var(wSamples)
	wmu = log( (wm^2)/(wv+wm^2)^0.5 )
	ws2 = log( 1 + (wv/(wm^2)) )
	#wupper1 = exp( wmu + (1.*(ws2^0.5)) )
	#wlower1 = exp( wmu - (1.*(ws2^0.5)) )
	#wupper2 = exp( wmu + (2.*(ws2^0.5)) )
	#wlower2 = exp( wmu - (2.*(ws2^0.5)) )
	wupper3 = exp( wmu + (3.*(ws2^0.5)) )
	wlower3 = exp( wmu - (3.*(ws2^0.5)) )



	#OUTPUT
	ding(1)
	layout( matrix(c(1, 2, 3), nrow=1, ncol=3) )	

	plot( seq(1, i+1), Zmax,
		"l",
		xlab="Sweeps",
		ylab="Value",
		main="Best Z Value"
	)
	abline( v=upper, col="red") 
	abline( v=lower, col="red") 
	abline( v=m, col='blue')

	L = length(seq(1, i))
	plot( seq(1, i), maxes[1:i], 
		main=sprintf('Sweeps: %d\n', i), 
		ylab='Max Improvement\n', 
		col=c(rep("black", L-1), "blue"),
		ylim=c(0, max(wupper3, maxes[1:i]))
	)	
        left = max((L-(W-1)), 1)
	#segments( left, wupper1, i, wupper1, col="green")
	#segments( left, wlower1, i, wlower1, col="green")
   	segments(    1, wupper3, left, wupper3, col="red", lty=2)
	segments(    1, wlower3, left, wlower3, col="red", lty=2)
	segments( left, wupper3,    i, wupper3, col="red")
	segments( left, wlower3,    i, wlower3, col="red")
	segments( left, wm, i, wm, col="black")
	
	hist( wSamples, add=F, freq=F,
		main="Window Samples"
	)
	abline( v=wupper3, col="red") 
	abline( v=wlower3, col="red")
	abline( v=m, col='blue' )
	abline( v=wm, col='black' )	
	legend("topright",
		legend=c("~99.7% Interval", "x-bar", "x-bar-bar"), 
		col=c("red", "blue", "black"),
		lty=1
	)

		
	#Sys.sleep(5)	
}



