rm(list=ls())

#
#FUNCTIONS
#

#
B = function(Bback, r, K, Cback){ 
	out = Bback + r*Bback*(1-(Bback/K)) - Cback 
	return(out)	
}

#
logIt = function(q, b){ log(q) + log(b) }

#
ll = function(p){
	#
	s  = p[1]
	q  = p[2]
	r  = p[3]
	K  = p[4]	
	#C0 = p[5]
	B0 = K #p[6]
	#
	Bnow = B(B0, r, K, C0) 
	out  = Y*log(2*pi*s^2)/2
	for(t in 1:Y){
		out  = out + ((log(y[t])-logIt(q, Bnow))^2)/(2*s^2)
		Bnow = B(Bnow, r, K, C[t])
	}
	#
	return(out)
}

#
#DATA
#

#
year = 1965:1987
y = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
C = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)

#
#MAIN
#

#misc
Y  = length(year)
C0 = 0#C[1]
#optimization
s0 = 0.1
q0 = 1
r0 = 1
K0 = 10000
#B0 = 1000
par = c(s0, q0, r0, K0)
optimOut = optim(par, ll)
#
sHat  = optimOut$par[1]
qHat  = optimOut$par[2]
rHat  = optimOut$par[3]
KHat  = optimOut$par[4]
B0Hat = KHat
#
Bnow = B(B0Hat, rHat, KHat, C0) 
IHat = matrix(NaN, nrow=1, ncol=Y)
BHat = matrix(NaN, nrow=1, ncol=Y+1)
colnames(IHat) = year
colnames(BHat) = c(year, year[Y]+1)
BHat[1] = Bnow
for(t in 1:Y){
	IHat[t]   = exp(logIt(qHat, BHat[t]))	
	BHat[t+1] = B(BHat[t], rHat, KHat, C[t])
}
#
plot(year, y)
lines(year, IHat)
#
writeLines(sprintf('rHat: %f', rHat))
writeLines(sprintf('KHat: %f', KHat))
writeLines(sprintf('qHat: %f', qHat))
writeLines(sprintf('sHat: %f', sHat))
writeLines(sprintf('B0Hat: %f', B0Hat))
