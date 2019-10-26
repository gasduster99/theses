library(tgp)

source('gifMake.r')
source('tSource.r')
source('alerts.r')

rosenbrock = function(x){
	out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
	return(out)
}

rosenbrockHD = function(x){
	d = dim(x)[2]
	out = 0
	for (i in seq(1, d-1)){
		inc = 100*(x[,i]^2 - x[,i+1])^2 + (x[,i] - 1)^2
		out = out + inc
	}
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

sortRes = function(res){
	M = max(res$rank)
	out = res
	for (i in seq(1,M)){ out[i,]=res[which(res$rank==i),] }	

	return(out)
}


#TUNING PARAMETERS

#initially I use 40/80
IN = 50
#sample size per iteration
n = 3
#sweeps of optimization
M = 60


#INITIALIZING
fName = 'rosenbrockHD' #'rosenbrock' #'rastrigin' #
f = eval(parse(text=fName)) #deparse(substitute(f))
#domain boundaries
rect = cbind(c(-5,-5,-5,-5), c(5,5,5,5)) #cbind(c(-3,-3), c(3,3)) #cbind(c(-2,-2), c(2, 2)) # cbind(c(-5,-5), c(5,5)) #cbind(c(-5,-5),c(5,5)) #cbind(c(-2, -2), c(2, 2)) #
m = dim(rect)[1]
name = sprintf('%sm%dM%dn%diN%d', fName, m, M, n, IN)#Big

#intially sampling the points 
X = lhs(IN, rect)
jay = dim(X)[2]
#making the "sample" by evaluating the function at the sample points
Z = f(X)

#ITERATING

#dev.new()
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

	#dev.off()	
	#dev.new(width=14, height=7)
	layout(matrix(c(1,2), nrow=1, ncol=2))
	plot(out$obj, layout='surf')
	plot(out$obj, as='improv', layout='as', main=sprintf('Expected Improvement\nSweep: %d\n', it), zlim=c(-0.001, 0.01))#add zlim

	return( list(X, Z, out, zOpt, improv) )
}
init = list(X, Z, out=NULL, zOpt=NULL, improv=numeric(0))

gifOut = gifMake(meat, init, seq(1, M), 1, name, width=1400)
out = gifOut[[3]]
zOpt = gifOut[[4]]
improv = gifOut[[5]]

#for(it in seq(1, M)){ init=meat(init, it) }
#out = init[[3]]
#zOpt = init[[4]]
#improv = init[[5]]

#w = 14
#h = 7
#dev.new(width=w, height=h)
#layout(matrix(c(1,2), ncol=2))
#plot(seq(1,M), zOpt, 'l', main=sprintf('Best Z %s', name), xlab='Sweep', ylab='z')
#plot(seq(1, M), improv, 'l', main=sprintf('Max Log Improv %s', name), xlab='Sweep', ylab='Max Log Improv')

ding()

bit = 1#readline('Save? (1=yes): ')
if ((as.numeric(bit)==1) | (bit=='yes') | (bit=='y')){
	#pdf(width=w, height=h)
	#layout(matrix(c(1,2), ncol=2))
	#plot(seq(1,M), zOpt, 'l', main=sprintf('Best Z %s', name), xlab='Sweep', ylab='z')
	#plot(seq(1, M), improv, 'l', main=sprintf('Max Log Improv %s', name), xlab='Sweep', ylab='Max Log Improv')
	#dev.off()

	save.image(sprintf('%s.RData', name))
}

#save.image(sprintf('%s.RData', name))
#ding()

















#for (it in seq(1, M)){
#	#the skafold section
#	out = optim.step.tgp(f, X=X, Z=Z, rect=rect, prev=out, verb=0)
#	X = rbind(X, out$X)
#	Z = c(Z, f(out$X))
#
#	progress = rbind(progress, out$progress)
#}

#writeLines( sprintf('\n____ Sweep: %d ____\n', it) )
	#print( sortRes(cbind(rFit$improv,XX)[rFit$improv$rank<=10,]) )
	#writeLines('\n')



	#my moving rect idea (move rect like pattern search)
	#centX = mean(X[,1])
	#centY = mean(X[,2])
	#rect = cbind(c(centX-1, centY-1), c(centX+1, centY+1))

	#XX = lhs(nn, rect)

	#the first bit in the paper
	#rFit = bgp(X, Z, XX, improv=c(1,n), verb=0)
	#X = cbind(rFit$improv,XX)[rFit$improv$rank<=n, seq(3, 3+(jay-1))]
	#Z = f(X)
