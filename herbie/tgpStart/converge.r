rm(list=ls())

library(tgp)

#interesting stuff from 35-40
#load('rastM40n3Small.RData')
#load('rastM80n3.RData')
#load('rastM80n1.RData')
#load('rastM100n1.RData')

#load('rastM100n1Big.RData')
#load('rastM200n1Big.RData')

#load('roseM100n1Big.RData')
#load('rosenbrockM80n1iN20.RData')
#load('rosenbrockHDm3M40n3iN50.RData')
#load('rosenbrockHDm4M40n3iN50.RData')
load('rosenbrockHDm4M60n3iN50.RData')

### TUNING PARAMETERS ###
#defines how precise to round stuff to and how to deal with numerical derivatives and stuff
prec = 6
#the sample size of of the null model
#en = 25
#a small variance
nullV = 10^-18

anova = function(m){
        #making the variable name a string variable
	inName = deparse(substitute(m))
        I=dim(m)[1]
        J=dim(m)[2]
	#
        X = matrix(0, nrow=I*J, ncol=J)
        Y = matrix(NaN, nrow=I*J, ncol=1)
        count=1
        for (j in seq(1, J)){
                for(i in seq(1, I)){
                        if ( ((count-1)%/%I+1)==j ){
                                X[count, j] = 1
                        }
			#
                        Y[count] = m[i, j]
                        count = count+1
                }
        }
	#
        sum = summary(aov(Y~X))
	return(sum)
}

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

	last = tail(breaks, 1)
	set = 0
	setCDF = 0

	nullA = 0.06
	nullB = 0.0005
	null = rgamma(length(set), nullA, scale=nullB)

	#I use length to get around numeric(0) quirks
	if (length(last)!=0){ 
		set = imp[seq(last, length(imp))]

		#aOut = anova(datM)
		#t.test
		#wilcox.test
	
		setCDF = ecdf(set)
			
		ksOut = ks.test(set, 'pgamma', nullA, nullB)
		print(ksOut)
	}
	
	#high values for "res" will make it it easier to identify small variances
	res = 1
	#print(vars^res)

	

	l=c(1,2,3,4,5,6,7,8)	
	#dev.off()
	#dev.new(width=w, height=h)
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


w=14
h=7.8
#dev.new(width=w, height=h)
#layout(matrix(c(1,2), ncol=2))
#plot(seq(1,M), zOpt, 'l', main=sprintf('Best Z %s', name), xlab='Sweep', ylab='z')
#plot(seq(1, M), improv, 'l', main=sprintf('Max Log Improv %s', name), xlab='Sweep', ylab='Max Log Improv')

dev.new(width=w, height=h)
for (m in seq(2, M)){
	conv(m)
	readline(sprintf('Sweep: %d\n', m))
}




#dev.new(width=14)
#layout(matrix(seq(1,20), ncol=5))
#for (i in seq(1, 20)){
#m = 20 #M)
#}
#load(datName)


