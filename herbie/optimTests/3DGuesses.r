rm(list=ls())

rosenbrock <- function(x){
        n <- length(x)
        out = sum(100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
        return(out)
}

rastrigin <- function(x){
        out = 10*length(x)+sum(x^2-10*cos(2*pi*x))
        return(out)
}


#h = 0.1
h = 0.01
xRose = seq(-5, 5, h)
yRose = seq(-50, 50, h)
ZRose = matrix( NaN, nrow=length(xRose), ncol=length(yRose) )
for (i in seq(1,length(xRose))){
        for (j in seq(1,length(yRose))){
                ZRose[i,j] = rosenbrock( c(xRose[i], yRose[j]) )
        }
}

xRast = seq(-5, 5, h)
yRast = seq(-5, 5, h)
ZRast = matrix( NaN, nrow=length(xRast), ncol=length(yRast) )
for (i in seq(1,length(xRast))){
        for (j in seq(1,length(yRast))){
                ZRast[i,j] = rastrigin( c(xRast[i], yRast[j]) )
	}
}


g0=c(-4,4)

pdf('roseGuesses.pdf')
filled.contour( xRose, yRose, ZRose,
                        color.palette = jet.colors,
                        xlab='X',
                        ylab='Y',
                        main='Rosenbrock',
                        #plot.axes={points(g0, pch=20, col='black')}
)
dev.off()

pdf('rastGuesses.pdf')
filled.contour( xRast, yRast, ZRast,
                        color.palette = jet.colors,
                        xlab='X',
                        ylab='Y',
                        main='Rastrigin',
                        #plot.axes={points(g0, pch=20, col='black')}
)
dev.off()
