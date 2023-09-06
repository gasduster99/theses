rm(list=ls())

library(latex2exp)

#
#FUNCTIONS
#

#
vb = function(x, kappa, max, x0=0){ max*(1-exp(-kappa*(x-x0))) }
pLaw = function(x, a, b){ a*x^b }

#
l2_dist = function(f, g, lower=-0.1, upper=10) {
	f_diff = function(x) (f(x) - g(x))^2
 	sqrt(integrate(f = f_diff, lower=lower, upper=upper, subdivisions=10^5)$value)
}

#
#MAIN
#

#
min = -0.7

#
age = seq(min, 10, 0.1)

#
Linf = 1
Lkap = 1
lengthVB = vb(age, Lkap, Linf, min)
a = 1
b = 3#
wtViaLen = pLaw(lengthVB, a, b)
#
weightVB = vb(age, Lkap*0.85, pLaw(Linf, a, b))

#
png("vbLW.png")
plot(age, lengthVB, type='l', col='blue', lty=3, lwd=3, ylim=c(0, pLaw(Linf, a, b)), ylab="Length or Weight", xlab="Age", main="Length/Weight Relationships")
lines(age, wtViaLen, col='black', lty=3, lwd=3)
#
lines(age, weightVB, lwd=3)
#
legend('bottomright', legend=c('VB Length', TeX('$a($VB Length$)^b$'), 'VB Weight'), col=c('blue', 'black', 'black'), lwd=3, lty=c(2, 2, 1))
dev.off()

#
#OPTIMIZE
#

#
par = c(Lkap, 0)
o = optim(par,
        function(x){
        #
        #a = x[1]
        #b = x[2]
        #c = x[3]
        #d = x[4]
        #
        l2_dist(
                function(y){
                        vb(y, x[1], pLaw(Linf, a, b), x0=x[2])
                },
                function(y){
			#
                        lengthVB = vb(y, Lkap, Linf, min)
			wtViaLen = pLaw(lengthVB, a, b)
                	#
			return(wtViaLen)
		}
        )
        }
)

#
png("vbOpt.png")
plot(age, lengthVB, type='l', col='white', lty=3, lwd=3, ylim=c(0, pLaw(Linf, a, b)), ylab="Weight", xlab="Age", main="Weight Relationships")
lines(age, wtViaLen, col='black', lty=3, lwd=3)
#
weightVBOpt = vb(age, o$par[1], pLaw(Linf, a, b), x0=o$par[2])
lines(age, weightVBOpt, lwd=3)
#
legend('bottomright', legend=c(TeX('$a($VB Length$)^b$'), 'VB Weight'), col=c('black', 'black'), lwd=3, lty=c(2, 1))
dev.off()


