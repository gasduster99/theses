rm(list=ls())

library(tgp)

source('gifMake.r')
source('tSource.r')
source('alerts.r')
source('funky.r')


#interesting stuff from 35-40
#load('rastM40n3Small.RData')
load('rastM80n3.RData')
#load('rastM80n1.RData')
#load('rastM100n1.RData')

#load('rastM100n1Big.RData')
#load('rastM200n1Big.RData')

#load('roseM100n1Big.RData')
#load('rosenbrockM80n1iN20.RData')
#load('rosenbrockHDm3M40n3iN50.RData')
#load('rosenbrockHDm4M40n3iN50.RData')
#load('rosenbrockHDm4M60n3iN50.RData')

### TUNING PARAMETERS ###
#defines how precise to round stuff to and how to deal with numerical derivatives and stuff
prec = 6
#the sample size of of the null model
#en = 25
#a small variance
nullV = 10^-18

getZones = function(sciVec, prec){
	#define the decision metric so as to always display in scientific notation
	options(scipen=-10)

	#use strings to truncate scientific notation numbers to the precision defined relative to the step size of numerical derivative
	truncVec = c(numeric(0))
	for (i in seq(1, length(sciVec))){	
		str = as.character(sciVec[i])
		nc = nchar(str)
		#first bit
		nastyBit = substr(str, 1, prec+2)
		if (substr(str, 1, 1)=="-"){ nastyBit=substr(str, 1, prec+3) }
		#last bit
		nastyBit = sprintf('%s%s', nastyBit, substr(str, nc-3, nc)) 
		
		truncVec = c(truncVec, as.numeric(nastyBit))
	}
	
	uTrunc = unique(truncVec)
	#indicies that belong to each segment (w/sep indicies)
	partI = list(numeric(0))
	#diff values that beloging to each segment (w/sep values)
	partD = list(numeric(0))
	for (u in seq(1, length(uTrunc))){
		partI[[u]] = which(truncVec==uTrunc[u]) 
		partD[[u]] = sciVec[partI[[u]]]
	}

	out = list(partI, partD)
	return(out)
}

conv = function(m){
	its = seq(1, m)
	imp = improv[its]
	lmOut = btlm(X=its, Z=imp, verb=0)	
	d_dxOut = d_dx.tgp(lmOut, its)
	dlmOutdX = d_dxOut[[1]]
	pOb = d_dxOut[[2]] 
	
	getOut = getZones(dlmOutdX, prec-1)
	I = getOut[[1]]
	vars = c(numeric(0))	
	breaks = c(numeric(0)) 
	for (ii in seq(1, length(I))){
		v = var(imp[I[[ii]]]) 
		vars =  c(vars, v)
		#keeping track of the break locations
		if (length(I[[ii]])==1 & its[I[[ii]]]!=max(its) & its[I[[ii]]]!=min(its)){ breaks=c(breaks, I[[ii]]) }
	}	

	#last = tail(breaks, 1)


	l=c(1,2,3,4,5,6,7,8)	
	options(scipen=0)
	layout(matrix(l, nrow=2, ncol=length(l)/2, byrow=T))
	
	#1
	plot(its, zOpt[its], 'l', ylab='z',main='Best Z')
	abline(h=0, col='red')
	abline(v=breaks, col='blue')
	#2
	plot(lmOut, center='km', layout='surf', main=sprintf('Sweeps: %d\n', m), ylab='Max Improvement\n')
	#3
	plot(its, dlmOutdX, ylim=c(0.5*-10^-2, 0.5*10^-2), ylab='Slope', main='1st Derivative')
	abline(h=0, col='red')
	abline(v=breaks, col='blue')
	#4
	plot(pOb, center='km', layout='as', as='ks2', ylim=c(0, 8*10^-5))
	abline(h=0, col='red')
	#5
	hist(set, xlim=c(0, 0.004), ylim=c(0, 2000), probability=T)
	curve(dgamma(x, nullA, scale=nullB), 0, 0.004, 1000, add=T, col='red')
	#6	
	plot(setCDF, xlim=c(0, 0.004))
	curve(pgamma(x, nullA, scale=nullB), 0, 0.004, 1000, add=T, col='red')
	#7
	curve( -sqrt(9-x^2)-2, -3, 3, 1000, ylim=c(-7, 5), xlim=c(-5, 5) )
	curve( sqrt( 1/4 -(x-2)^2 ) + 2, -5, 5, add=T )
	curve( -sqrt(1/4 - (x-2)^2 ) + 2, -5, 5, add=T )
	curve( sqrt( 1/4 -(x+2)^2) + 2, -5, 5, add=T )
	curve( -sqrt( 1/4 -(x+2)^2 ) + 2, -5, 5, add=T )
	#8
	curve( -sqrt(9-x^2)-2, -3, 3, 1000, ylim=c(-7, 5), xlim=c(-5, 5) )
	curve( sqrt( 1/4 -(x-2)^2 ) + 2, -5, 5, add=T )
	curve( -sqrt(1/4 - (x-2)^2 ) + 2, -5, 5, add=T )
	curve( sqrt( 1/4 -(x+2)^2) + 2, -5, 5, add=T )
	curve( -sqrt( 1/4 -(x+2)^2 ) + 2, -5, 5, add=T )
}

d_dx = function(fun, x0, h=10^(-prec), order=4){ 
        #fun is a function to take the numerical derivative 
        #x0 is the stop on the domain where you should evaulate you derivative
        #h is the step size to use to take the numerical derivative (should be small, but 10^-6 is about as small as I can make it before numerical issues.

	#first order
        if (order==1){ out=(fun(x0+h)-fun(x0))/h}
	
	#fourth order
	if (order==4){
		out = (fun(x0-2*h)-fun(x0+2*h)) - 8*(fun(x0-h)-fun(x0+h))
		out = out/(12*h)
	}

        return(out)
}

d_dx.tgp = function(tgpOb, x0, h=10^(-6), order=4){
	#tgpOb is an tgp object outputted by a b* function 
        #x0 is the stop on the domain where you should evaulate you derivative
        #h is the step size to use to take the numerical derivative (should be small, but 10^-6 is about as small as I can make it before numerical issues.
	#order is the order of the numerical derivative used (1 or 4)

	pX = predict(tgpOb, XX=x0)

	#first order
        if (order==1){ 	
		pXH = predict(tgpOb, XX=x0+h)
		out = ( pXH$ZZ.km-pX$ZZ.km )/h
	}
	
	#fourth order
	if (order==4){
		pXh = predict(tgpOb, XX=x0-h)
		pXH = predict(tgpOb, XX=x0+h)
		pXHH = predict(tgpOb, XX=x0+2*h)
		pXhh = predict(tgpOb, XX=x0-2*h)		

		out = (pXhh$ZZ.km-pXHH$ZZ.km) - 8*(pXh$ZZ.km-pXH$ZZ.km)
		out = out/(12*h)
	}

        return(list(out, pX))
}

#TUNING PARAMETERS

#initially I use 40/80
IN = 50
#sample size per iteration
n = 3
#predictive sample size
nn = 200
#sweeps of optimization
M = 60


#INITIALIZING
fName = 'rosenbrockHD' #'rosenbrock' #'rastrigin' #
f = eval(parse(text=fName)) #deparse(substitute(f))
#domain boundaries
rect = cbind(c(-5,-5,-5,-5), c(5,5,5,5)) #cbind(c(-3,-3), c(3,3)) #cbind(c(-2,-2), c
m = dim(rect)[1]
name = sprintf('%sm%dM%dn%diN%d', fName, m, M, n, IN)#Big

#intially sampling the points 
X = lhs(IN, rect)
jay = dim(X)[2]
#making the "sample" by evaluating the function at the sample points
Z = f(X)

#out$obj$improv[out$obj$improv$rank==1,]
sortRes = function(res){
        M = max(res$rank)
        out1 = res
	print(M)
        for (i in seq(1,M)){ out1[i,]=res[res$rank==i,] }
	out1 = out1[seq(1, M),]
	
        return(out1)
}

meat = function(inits, it){
        X = inits[[1]]
        Z = inits[[2]]
        out = inits[[3]]
        zOpt = inits[[4]]
        improv = inits[[5]]

        out = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=out, verb=0, improv=c(1,n))
        ex = matrix(out$X, ncol=m)
        fex = f(ex)
        X = rbind(X, ex)
        Z = c(Z, fex)
        prog = out$progress
        improv = c(improv, prog$improv)	
        zOpt = c(zOpt, min(Z))
	
	I = out$obj$improv	
	EI = I$improv
	print( length(EI) )	

	#M = max(EI$rank)
	#for (i in seq(1, M)){ print(EI[EI$rank==i,]) }

	#EI = sortRes(out$obj$improv)
	#print(EI)

        return( list(X, Z, out, zOpt, improv) )
}
init = list(X, Z, out=NULL, zOpt=NULL, improv=numeric(0))

for(it in seq(1, M)){ init=meat(init, it) }
out = init[[3]]
Opt = init[[4]]
improv = init[[5]]


w=14
h=7.8

dev.new(width=w, height=h)
for (m in seq(2, M)){
	conv(m)
	readline(sprintf('Sweep: %d\n\n', m))
}















