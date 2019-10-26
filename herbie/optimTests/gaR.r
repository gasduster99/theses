library("GA")
#ga finds maximum values so to find minima you have to look at the negative 

f = function(x){
	out = sin(x) + 0.001*x
	return(out)
}

monitor = function(obj){
	curve(f, min, max, main=obj@iter, font.main = 1)
	points(obj@population, obj@fitness, pch = 20, col = 2)
	rug(obj@population, col = 2)
	#Sys.sleep(0.2)
}

min = 0
max = 4*pi

pdf('gaSearchingR.pdf', width=14, height=9)
layout(matrix(seq(1, 50), 5, 10, byrow=T))
GA = ga(type = "real-valued", fitness = f, min = min, max = max, maxiter=50, monitor=monitor)
dev.off()

print( summary(GA) )

pdf('gaConvergeR.pdf')
plot(GA)
dev.off()











#layout( matrix(c(1,2), 1, 2, byrow=TRUE) )
