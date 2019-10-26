rastrigin <- function(x){
	out = 10*length(x)+sum(x^2-10*cos(2*pi*x))
	return(out)
}

rosenbrock <- function(x){
	n <- length(x)
	out = sum(100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
	return(out)
}


h=0.1
xDom = seq(-5, 5, h)
yDom = seq(-5, 5, h)
Zrasta = matrix( NaN, nrow=length(xDom), ncol=length(yDom) )
Zrose = matrix( NaN, nrow=length(xDom), ncol=length(yDom) )

for (i in seq(1,length(xDom))){
	for (j in seq(1,length(yDom))){
		Zrasta[i,j] = rastrigin( c(xDom[i], yDom[j]) )
		Zrose[i,j] = rosenbrock( c(xDom[i], yDom[j]) )
	}
}

persp(xDom, yDom, Zrasta)
dev.new()
persp(xDom, yDom, Zrose)

