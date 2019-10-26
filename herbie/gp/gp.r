rm(list=ls())

library(MASS)

K = function(a, b){
	#Give this function 2 vectors a and b
	#a sits on the rows
	#b sits on the columns
	out = matrix(NaN, nrow=length(a), ncol=length(b))
	for (i in seq(1, length(a))){
		for (j in seq(1, length(b))){
			d = a[i]-b[j]
			out[i,j] = exp(-(d/4)^2)
		}
	}
	
	return(out)	
}
Knug = function(a, b, nug){ 
        #Give this function 2 vectors a and b
        #a sits on the rows
        #b sits on the columns

        out = matrix(NaN, nrow=length(a), ncol=length(b))
        for (i in seq(1, length(a))){
                for (j in seq(1, length(b))){
                        d = a[i]-b[j]
                        out[i,j] = exp(-(d/4)^2) + nug
                }   
        }   
     
        return(out)     
}


#layout(matrix(c(1,2), nrow=1, ncol=2))

fullData = read.table('data.txt')
x = fullData$V1
z = fullData$V2
#
l = length(fullData$V2)
p = 5 
keep = rep(1, l)
keep[seq(1, l, p)] = 0
sd = 0.015
#
z0 = fullData$V2[seq(1, l, p)]#+rnorm(l, 0, sd)#*keep
z1 = fullData$V2[seq(1, l, p)]#+rnorm(l, 0, sd)#*keep
z2 = fullData$V2[seq(1, l, p)]#+rnorm(l, 0, sd)#*keep
#
SIG0 = K(x[seq(1, l, p)],x[seq(1, l, p)])
SIG = K(x,x)
SIGI0 = solve(SIG0)
SIGI = solve(SIG)
#
domain = seq(0, 25, 0.01)
so = K(x, domain)
so0 = K(x[seq(1, l, p)], domain)
y = t(z)%*%SIGI%*%so
#
y0 = t(z0)%*%SIGI0%*%so0
y1 = t(z1)%*%SIGI0%*%so0
y2 = t(z2)%*%SIGI0%*%so0
#plot(fullData$V1, fullData$V2) #z)#
pdf('surrogateStep.pdf')
plot(domain, y, 'l', 
	lwd = 5,
	xlim = c(fullData$V1[2], fullData$V1[l-1]),
	ylim = c(-2, 3),
	xlab='x',
	main='One Surrogate Modeling Iteration'
	)
#lines(domain, y0, lty=2)
#lines(domain, y1, lty=2)
lines(domain, y2, col='red', lwd=2, lty=2)
points(seq(1, l, p), fullData$V2[seq(1, l, p)], pch=19, col='red')
segments(domain[y2==min(y2)], -3, domain[y2==min(y2)], min(y2), 0.05, lwd=2, col='red', lty=2)
segments(domain[y==min(y)], -3, domain[y==min(y)], min(y), lwd=5)
legend('bottomright', c('Truth', 'Surrogate Model'), 
	lty=c(1, 2),
	lwd=c(5,2), 
	col=c('black', 'red')
)
dev.off()
#color = colorRampPalette(c("orange", "red", "violet", "blue", "lightblue", "lightgreen", "darkgreen")) 
#try = seq(0, 0.2, 0.001)
#for (n in try){
#	SIG = Knug(x, x, n)
#	SIGI = solve(SIG)
#	
#	domain = seq(0, 25, 0.01)
#	so = Knug(x, domain, n)
#	y = t(z)%*%SIGI%*%so
#	
#	#plot(fullData$V1, fullData$V2)
#	lines(domain, y, col=color(length(try)))
#}
