rm(list=ls())

#
library(Ryacas)
library(pracma)
library(stringr)

#
#FUNCTIONS
#

#prepmath
##
#R = ysym("alpha*B*(1-gamma*beta*B)^(1/gamma)")
##
#Nbar = ysym("R - (M+F)*N")
#Nbar = with_value(Nbar, "R", R)
#Nbar = ysym("(alpha*B*(1-gamma*beta*B)^(1/gamma))/(M+F)") #solve(Nbar, 'N')[1]
##
##Bbar = ysym("w*R + k*(W*N-B) - (M+F)*B")
##Bbar = with_value(Bbar, "R", R)
##Bbar = with_value(Bbar, "N", Nbar)

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%28W*a*B*%281-gamma*bet
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
BBar = function(FF, M, k, w, W, alpha, beta, gamma){ eval(as_r( with_value(Bbar, "F", ysym("FF")) )) }
#
FBbar = ysym("F*B")
FBbar = with_value(FBbar, "B", Bbar)
FBbar = with_value(FBbar, "F", ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
FMsy = function(M, k, w, W, alpha, beta, gamma){ 
	uniroot(function(FF){ eval(as_r(dFBdF)) }, c(0, 10))$root 
}
#beta does not matter for either of these
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma){
	#FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(alpha){ eval(as_r(dFBdF)) }, c(eps(), 100))$root
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta){
	#FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(as_r(dFBdF)) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
	f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
	uniroot(f, c(0, 10), tol=eps())$root
}

#Does BStar/BZero depend on both alpha and gamma 

#I think this is not th ebest way to get alpha 
#getAlphaZeta = function(zeta, FF, M, k, w, W, beta, gamma){
#	f = function(a){ 
#		BBar(FF, M, k, w, W, a, beta, gamma)/BBar(0, M, k, w, W, a, beta, gamma) - zeta
#	}
#	uniroot(f, c(eps(), 10))$root
#}
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
	f = function(g){ 
		BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - zeta
	}
	uniroot(f, c(-10, 10))$root
}

#
getZetaBH = function(x, M, k, W, a0){
        #
	w = W*(1-exp(-k*a0))
	#
        gamma = -1
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaBH = Vectorize(getZetaBH, "x")

#
#SYMBOL
#

#
w = function(a, wi, k){ wi*(1-exp(-k*a)) }
#
M = 0.2
k = 0.2
#
a = 2
W = 1
w = w(a, W, k)
#
alpha = 5
beta = 1
gamma = -0.5
#
Fmsy = FMsy(M, k, w, W, alpha, beta, gamma) 
xi = Fmsy/M

#
BStar = BBar(Fmsy, M, k, w, W, alpha, beta, gamma) 
BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
zeta = BStar/BZero

#
xi = 1
zeta = 1/(xi+2)
B0 = 10000
#

#AN IMPLICATE IDEA
xTol = 10^-3
zTol = xTol
bTol = 10
xiHat = xi+xTol*10
zetaHat = zeta+zTol*10
#
xis = c()
zetas = c()
zeros = c()
#for(i in 1:10){
i = 0
#| abs(BZero-B0)>=bTol
while( abs(xiHat-xi)>=xTol | abs(zetaHat-zeta)>=zTol ){
	#
	gamma = getGammaZeta(zeta, xi*M, M, k, w, W, alpha, beta)
	alpha = getAlphaFmsy(xi*M, M, k, w, W, beta, gamma)
	beta  = getBeta(B0, M, k, w, W, alpha, gamma)
	#
	BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
	xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
	zetaHat = BBar(xiHat*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
	#
	xis = c(xis, xiHat)
	zetas = c(zetas, zetaHat)
	#
	i = i+1
}


#BH RPS
f = function(x){1/(x+2)}
curve(f(x), 0.01, 3, lwd=3)
curve(getZetaBH(x, M, k, W, a), 0.01, 3, add=T)
curve(getZetaBH(x, M, 1, W, a), 0.01, 3, col=2, add=T)
curve(getZetaBH(x, M, 0.15, W, a), 0.01, 3, col=3, add=T)

##the constraint is invariant to W
#dev.new()
#curve(getZetaBH(x, M, k, W, a), 0.01, 3)
#curve(getZetaBH(x, M, k, 5, a), 0.01, 3, col=2, add=T)
#curve(getZetaBH(x, M, k, 0.5, a), 0.01, 3, col=3, add=T)

#
dev.new()
f = function(x){1/(x+2)}
curve(f(x), 0.01, 3, lwd=3)
curve(getZetaBH(x, M, k, W, a), 0.01, 3, add=T)
curve(getZetaBH(x, M, k, W, 10), 0.01, 3, col=2, add=T)
#curve(getZetaBH(x, M, k, W, a), 0.01, 3, col=3, add=T)



##
#zetaBH = c()
#xiBH = seq(0.1, 4, 0.1)
#for(x in xiBH){
#	#
#	gamma = -1 #getGammaZeta(zeta, xi*M, M, k, w, W, alpha, beta)
#	alpha = getAlphaFmsy(x*M, M, k, w, W, beta, gamma)
#	beta  = getBeta(B0, M, k, w, W, alpha, gamma)
#	#
#	BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
#	xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
#	zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
#	#
#	zetaBH = c(zetaBH, zetaHat)
#}

##
#BBH = with_value(Bbar, "gamma", -1)
##
#FBBH = ysym("F*B")
#FBBH = with_value(FBBH, "B", BBH)
#FBBH = with_value(FBBH, "F", ysym("FF"))
##
#dFBdFBH = deriv(FBBH, "FF")
#yac('MaxEvalDepth(100000)')
#FBH = solve(dFBdFBH, 'FF')



##dig into LHS code
#	#I think sample alpha with a T
#	#get gamma
#	#get beta
#	#get Fmsy
#	#figure out which zeta was sampled
#	#
#
##
#f = function(alpha, gamma, x, z){
#	#
#	gamma = getGammaZeta(zeta, xi*M, M, k, w, W, alpha, beta)
#	alpha = getAlphaFmsy(xi*M, M, k, w, W, beta, gamma)
#	beta  = getBeta(B0, M, k, w, W, alpha, gamma)
#	#
#	xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
#	zetaHat = BBar(xiHat*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
#	#
#	return( (xiHat-x)^2 + (zetaHat-z)^2 )
#}
##
#optOut = optim(c(alpha, gamma), function(a, b){f(a, b, xi, zeta)}, 
#	lower = c(eps(), -10),
#	upper = c(10, 10),
#	method = "L-BFGS-B"
#)






##
#R = function(B, alpha, beta, gamma){
#	alpha*B*(1-gamma*beta*B)^(1/gamma)
#}

#Fmsy = uniroot(function(x){dYdF(x, M, k, w, W, alpha, beta, gamma)}, c(0, 10))$root

##yac('MaxEvalDepth(100000)')
##Fmsy  = solve(dFBdF, "F")

#
#yac_expr(dFBdF)


##
#R = ysym("alpha*B*(1-gamma*beta*B)^(1/gamma)")
##
#Nbar = ysym("R - (M+F)*N")
#Nbar = with_value(Nbar, "R", R)
#Nbar = ysym("(alpha*B*(1-gamma*beta*B)^(1/gamma))/(M+F)") #solve(Nbar, 'N')[1]
##
#Bbar = ysym("w*R + k*(W*N-B) - (M+F)*B")
#Bbar = with_value(Bbar, "R", R)
#Bbar = with_value(Bbar, "N", Nbar)

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%28W*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%29%2F%28M%2BF%29-B%29-%28M%2BF%29*B+for+B&assumption=%7B%22C%22%2C+%22gamma%22%7D+-%3E+%7B%22Variable%22%7D
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
#Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
##
#FBbar = ysym("F*B")
#FBbar = with_value(FBbar, "B", Bbar)
#FBbar = with_value(FBbar, "F", ysym("FF"))
##
#dFBdF = deriv(FBbar, "FF")
##yac('MaxEvalDepth(10000)')
##Fmsy = solve(dFBdF, 'FF')

##
#dYdF = function(FF, M, k, w, W, alpha, beta, gamma){ eval(as_r(dFBdF)) }




##
#Bbar = ysym("R * (w+k*W/(M+F)) / (M+F+k)")
#Bbar = with_value(Bbar, "R", R)
##NOTE: need to solve above for 'B'
#



##
#R = ysym("alpha*P/(1+beta*P^(1/gamma))")
#dPdt = ysym("R - (M+F)*P")
#dPdt = with_value(dPdt, "R", R)
##
#solve(dPdt, "P")

##
#expr = yacas("alpha*P/(1+beta*P^(1/gamma))-(M+F)*P==0")
#sol = Solve(expr, "P")

#yac_expr("MaxEvalDept(1000000)")
#dPdt = ysym("alpha*P/(1+beta*P^(1/gamma))-(M+F)*P")
#Pbar = solve(dPdt, "P")#, MaxEvalDept=10^4)

##
#P = ysym("(alpha/(M+F)-1)^gamma * beta^-gamma")
#L = ysym("alpha*(M+F)-(M+F)^2-alpha*gamma*F")
##
#alpha=1
#beta=1
#gamma=1
#M=0.2
##
#FStar = solve(L, "F")
#st = deparse(as_r(FStar))
#rm = str_replace( st, "F == ", "")
#l = eval(eval(parse(text=rm)))
#rep = str_replace(deparse(as_r(FStar[l>0])), "F == ", "")
#rep = str_replace( rep, "sqrt", "Sqrt")
#fStar = ysym(eval(parse(text=rep)))
##f = ysym("-(alpha*gamma-(alpha-2*M)-Sqrt((alpha-2*M-alpha*gamma)^2+4*(alpha*M-M^2)))/2")
#
###FStar = as_y(rm[l>0])
#
##
#PStar = with_value(P, "F", fStar)
#P0 = with_value(P, "F", 0)
##
#M = c(fStar, PStar, P0)
#abg = solve(M, c(0.2, 2000, 3000), c("alpha", "beta", "gamma"))
#
##M = with_value(M, "M", 0.2)
##abg = solve(M, c(0.2, 2000, 3000), c("alp", "bet", "gam"))


#
#JUNK
#

#deparse(as_r(FStarMABG[2]))

#FStarM = with_value(FStar, "M", 0.2)
#FStarMA = with_value(FStarM, "alpha", 1)
#FStarMAB = with_value(FStarMA, "beta", 1)
#FStarMABG = with_value(FStarMAB, "gamma", 1)

#dPdF = deriv(P, "F")
#dYdF = with_value(ysym("P+F*dPdF"), "P", P)
#dYdF = with_value(dYdF, "dPdF", dPdF)

