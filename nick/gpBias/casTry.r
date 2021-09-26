rm(list=ls())

#
library(Ryacas)
library(stringr)

#
#
#

#

#
#SYMBOL
#

#
R = ysym("alpha*P/(1+beta*P^(1/gamma))")
dPdt = ysym("R - (M+F)*P")
dPdt = with_value(dPdt, "R", R)
#
solve(dPdt, "P")

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

