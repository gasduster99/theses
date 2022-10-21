rm(list=ls())

#
library(Ryacas)
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

