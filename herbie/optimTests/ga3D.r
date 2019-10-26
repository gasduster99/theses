library("GA")

alpha=c('aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al','am','an','ao','ap','aq','ar','as','at','au','av','aw','ax','ay','az','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z')

rastrigin <- function(x){
        out = 10*length(x)+sum(x^2-10*cos(2*pi*x))
        return(out)
}

rosenbrock <- function(x){
        n <- length(x)
        out = sum(100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
        return(out)
}

rastrigin <- function(x){
        out = 10*length(x)+sum(x^2-10*cos(2*pi*x))
        return(out)
}

f = function(x){
        out = sin(x) + 0.001*x
        return(out)
}

monitor = function(obj){
	jpeg(sprintf('gaSearchingRasta/%s.jpeg', alpha[obj@iter]))#, width=10, height=90/14)
	#dev.new(width=14, height=9)
	#layout(matrix(seq(1, 5), 5, 10, byrow=T))
	#par(mar=c(1,1,1,1))
	filled.contour( xCon, yCon, ZCon,
                	color.palette = jet.colors,
                	xlab='X',
                	ylab='Y',
                	main=obj@iter,
                	plot.axes={points(obj@population, pch=20, col='black')}
	)
	dev.off()
        #curve(f, min, max, main=obj@iter, font.main = 1)
        #points(obj@population, obj@fitness, pch = 20, col = 2)
        #rug(obj@population, col = 2)
        #Sys.sleep(0.2)
}


#h = 0.1
h = 0.01
xPer = seq(-5, 5, h) #seq(-2, 2, h)
yPer = seq(-5, 5, h) #seq(-3, 5, h)
ZPer = matrix( NaN, nrow=length(xPer), ncol=length(yPer) )
for (i in seq(1,length(xPer))){
        for (j in seq(1,length(yPer))){
                ZPer[i,j] = rastrigin( c(xPer[i], yPer[j]) ) #rosenbrock( c(xPer[i], yPer[j]) )
        }
}

h = 0.1
xCon = seq(-5, 5, h)
yCon = seq(-5, 5, h)
ZCon = matrix( NaN, nrow=length(xCon), ncol=length(yCon) )
for (i in seq(1,length(xCon))){
        for (j in seq(1,length(yCon))){
                ZCon[i,j] = rastrigin( c(xCon[i], yCon[j]) )  #rosenbrock( c(xCon[i], yCon[j]) )
        }
}

#persp3D(xDom, yDom, Zrasta)
pdf('rastPersp.pdf')
persp3D(xPer, yPer, ZPer,
	theta=100,
	phi=30,
	xlab='X',
	ylab='Y',
	zlab='Z', 
	main='Rastrigin Function'
	)
dev.off()

pdf('rastContour.pdf')
filled.contour( xPer, yPer, ZPer, color.palette=jet.colors, xlab='X', ylab='Y', main='Rastrigin Function' )
dev.off()

#min = c(-5, -5)
#max = c(5, 5)
#GA = ga( type = "real-valued",
#	 fitness = function(x){ -rastrigin(c(x[1], x[2])) },
#	 min = min, max = max,
#	 maxiter=50,
#	 #elitism=,
#	 monitor=monitor
#	)
#
#print( summary(GA) )
#
#pdf('gaConverge3DRasta.pdf')
#plot( GA )
#dev.off()








