#rm(list=ls())

#
#FUNCTIONS
#

#BH
bev = function(f){ f/(f*(f+2)) }

#
disX = function(x, p){
	#x: a value on the x-axis
	#p: a vector of reference point
	((p[1]-x)^2+(p[2]-bev(x))^2)^0.5
}

#
missSpec = function(x, y, xlim){
	#x: x value of 'truth' point
	#y: y value of 'truth' point
	#xlim: limits of search domain
	optimise(disX, xlim, p=c(x, y)
}

#
pCircle = function(c, r){
	#c: center
	#r: radius
	xs = seq(c[1]-r, c[1]+r, 0.00001)#0.01)#	
	#r^2 = x^2+y^2
	yp = (r^2-(xs-c[1])^2)^0.5
	ym = -yp
	#plot
	lines(xs, yp+c[2])
	lines(xs, ym+c[2])
	segments(xs[1], ym[1]+c[2], xs[1], ym[1]+c[2])
	segments(xs[length(xs)], ym[length(ym)]+c[2], xs[length(xs)], yp[length(yp)]+c[2])
}


##
##TEST
##
#
##
#curve(bev, 0, 5, ylim=c(0, 1), col='red', lwd=3)
##
#xs = seq(0.5, 3, 0.5)
#ys = seq(0.2, 0.6, 0.1)
#for(x in xs){
#	for(y in ys){
#		#
#		op  = optimise(disX, c(0, 5), p=c(x, y))
#		#
#		points(x, y, pch=19, col='blue')
#		pCircle(c(x, y), op[[2]])
#		segments(x, y, op[[1]], bev(op[[1]]), col='blue')
#	}
#}











##BH function
#dash = function(t){ t/(t*(t+2)) }
##distamce from (x, dash(x)) to (1, 0.6) #sqrt( (x-dot[1])^2 + (dash(x)-dot[2])^2 ) }
##dist = function(x){ norm( cbind(dot, c(x, dash(x))), "2" ) }
##dist = function(x, p){ norm( cbind(p, p-c(x, dash(x))), "2" ) }
#dist = function(x){ sqrt(sum((c(x, dash(x))-p)^2)) }
#
##
#curve(dash, 0, 5, ylim=c(0, 1))
##
#xs = seq(0, 3, 0.5)
#ys = seq(0.2, 0.6, 0.1)
#for(x in xs[1]){
#	for(y in ys[1]){
#		#point of interest
#		p = c(x, y)
#		#op  = optimize(dist, c(0, 5))
#		op = optimize(dist, c(0, 5))
#		#
#		segments(x, y, op[1][[1]], dash(op[1][[1]]))
#	}
#}
