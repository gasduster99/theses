rm(list=ls())

#
library(plgp)
library(crch)
library(gMOIP)
library(Ryacas)
library(pracma)
library(stringr)
library(rootSolve)
library(latex2exp)
library(RColorBrewer)

#
source('ddClass0.1.1.r')

#
#DELAY FUNCTIONS
#

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
Bbar_r = as_r( with_value(Bbar, "F", ysym("FF")) )
BBar = function(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#
FBbar = ysym("F*B")
FBbar = with_value(FBbar, "B", Bbar)
FBbar = with_value(FBbar, "F", ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
dFBdF_r = as_r(dFBdF)
FMsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10), tol=.Machine$double.eps)$root
}
#
NBar = function(FF, M, k, w, W, alpha, beta, gamma){
        #
        BB = BBar(FF, M, kappa, ww, WW, alpha, beta, gamma)
        (alpha*BB*( 1-beta*gamma*BB )^(1/gamma))/(M+FF)
}
NCar = function(x, FF, M, k, w, W, alpha, beta, gamma){
        #
        #BB = BBar(FF, M, kappa, ww, WW, alpha, beta, gamma)
        (alpha*x*( 1-beta*gamma*x )^(1/gamma))/(M+FF)
}
#beta does not matter for either of these
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma, dFBdFexpr=dFBdF_r){
        #
        #capture.output(
        out <- tryCatch({
                uniroot(function(alpha){ eval(dFBdFexpr) }, c(eps(), 100))$root
        }, error=function(err){
                out = NA
        })#, file="/dev/null")
        #
        return(out)
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta, dFBdFexpr=dFBdF_r){
        #FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(dFBdFexpr) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
        f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
        uniroot(f, c(0, 10), tol=eps())$root
}

#beta determines Bzero
getBetaBmsy = function(BB, M, k, w, W, alpha, gamma){
        f = function(b){ 
		FF = FMsy(M, k, w, W, alpha, b, gamma)
		BBar(FF, M, k, w, W, alpha, b, gamma) - BB
	}
        uniroot(f, c(0, 10), tol=eps())$root
}

#gamma|alpha, zeta
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
        f = function(g){
                BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - ze
        }
        uniroot(f, c(-10, 10))$root
}
#
getZeta = function(FF, M, k, w, W, alpha, beta, gamma){
        BBar(FF, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
}

#
invert = function(zeta, xi, B0, M, k, w, W, alphaStart=1, betaStart=1, gammaStart=1){
        #
        alpha = alphaStart
        beta  = betaStart
        gamma = gammaStart
        #
        xTol = 10^-3
        zTol = xTol
        bTol = 10
        xiHat = xi+xTol*10
        zetaHat = zeta+zTol*10
        #
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
                #xis = c(xis, xiHat)
                #zetas = c(zetas, zetaHat)
                #i = i+1
        }
        #
        return( data.frame(alpha=alpha, beta=beta, gamma=gamma) )
}

#
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#
getZetaBH = function(x, M, k, W, aS, a0){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
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
getZetaG = function(x, M, k, W, aS, a0, G){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
        #
        gamma = G
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaG = Vectorize(getZetaG, "x")

#DIFFERENTIAL EQUATION 

#
der = function(t, Y, lalpha, lbeta, gamma, aS, a0, WW, kappa, catch, B0){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        N = Y[1]
        B = Y[2]
        #
        if( (t-aS)<1){ #0.1
                Blag = B0
        }else{
                Blag = lagvalue(t-aS)[2]
        }
        #
        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        alpha = exp(lalpha)
        beta = exp(lbeta)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
        #
        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
        FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
        #
        out = c(N=NA, B=NA)
        out[1] = R - (M+FF)*N
        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B
        #
        return( list(out, ww*R, kappa*(WW*N-B)) )
}

#
SRR = function(B, mod){
        #
        alpha = exp(mod$lalpha)
        beta = exp(mod$lbeta)
        gamma = mod$gamma
        #
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
        #
        return(R)
}
SRR = Vectorize(SRR, 'B')

#
surplus = function(B, mod){
        #
        alpha = exp(mod$lalpha)
        beta = exp(mod$lbeta)
        gamma = mod$gamma
        #
        ww = vbGrow(mod$aS, mod$kappa, mod$WW, mod$a0)
        #
        FF = 0 #dat$FMsy
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
        N = R/(mod$M+FF)
        #Y = ww*R + mod$kappa*(mod$WW*N-B) - (M+FF)*B
        #Y = ww*R + mod$kappa*(mod$WW*R/(M+FF)-B) - (M+FF)*B
        Y = ww*R + mod$kappa*mod$WW*N - mod$kappa*B - mod$M*B - FF*B
        #fDat = function(x){ wwDat*SRR(x, dat) -dat$M*x-dat$kappa*x + 
        #       dat$kappa*dat$WW*NCar(x, 0, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma) }
        # (alpha*x*( 1-beta*gamma*x )^(1/gamma))/(M+FF)
        return(Y)
}
surplus = Vectorize(surplus, 'B')

#CATCH

#
fContrast = function(con){
        TT = 31 #length(hake)
        tt = 1
        time = tt:TT
        #
        mid = round(TT/2)
        cMax = 1.6^con #2^con  
        cMin = 0.4^con #0.3^con #0.2^con
        b = log(cMin/cMax)/(tt-mid)
        a = exp(log(cMax)-b*mid)
        rSlope = (cMax-1)/(mid-1)
        rb = 2*rSlope*mid + cMax-rSlope*mid
        FtFmsy = (a*exp(b*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
        #
        FtFmsy = c(FtFmsy, rep(1, 15))
}

#
#MAIN
#

#
path = "./modsDDExpT45N150A0-1AS2K0.1/"
startWho = sprintf("%s/fit3KA_xi1.463_zeta0.418.rda", path)
whoMight = as.data.frame(t(list.files(path=path, pattern=glob2rx("fit3KA*.rda"))))
colnames(whoMight)=whoMight
dat = readRDS(startWho)
dat$whoAmI = startWho
#aS = dat$aS
B0 = 10000
WW = 1
M = dat$M
dat$xi = NULL
dat$zeta = NULL
#
dat$catch = fContrast(0)
#
cols = brewer.pal(9, 'Set1')

#
png('growthTriptic.png', width=480*3)
layout(t(1:3))

#
#BH
#

#Growth
dat$aS = 1 #4 #10
dat$kappa = 0.5 #0.2 #0.1
#Recruitment
dat$alpha = 1.2
dat$lalpha= log(dat$alpha)
#dat$beta  = 0.0002
#dat$lbeta = log(dat$beta)
dat$gamma = -1
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.00006 #-0.00007
dat$lbeta = log(dat$beta)
dat$iterate()
BB = dat$B[46]
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Beverton-Holt", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=2, col=cols[3], cex.main=3, cex.axis=2, cex.lab=3, lty=2)
#toned down
dat$aS = 2
dat$kappa = 0.1
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
dat$lbeta = log(dat$beta)
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[3], lty=3)
#prod
dat$aS = 0.1
dat$kappa = 10
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
dat$lbeta = log(dat$beta)
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[3], lty=1)
#legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[3], lwd=2, lty=rev(1:3), cex=2)
legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=1$,    $\\kappa=0.5$, $w(a_s)\\approx %s$", round(vbGrow(1, 0.5, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[3], lwd=2, lty=rev(1:3), cex=2)

#
#RICKER
#

#Growth
dat$aS = 1 #4 #10
dat$kappa = 0.5 #0.2 #0.1
#Recruitment
dat$alpha = 1.2 
dat$lalpha= log(dat$alpha)
#dat$beta  = 0.0002 
#dat$lbeta = log(dat$beta)
dat$gamma = 0.0001
dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)#-0.000015
dat$lbeta = log(dat$beta)
dat$iterate()
BB = dat$B[46]
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Ricker", ylim=c(0,B0), xlab="Time", ylab="Biomass", col=cols[2], lwd=2, cex.main=3, cex.axis=2, cex.lab=3, lty=2)
#toned down
dat$aS = 2
dat$kappa = 0.1
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
dat$lbeta = log(dat$beta)
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[2], lty=3)
#prod
dat$aS = 0.1
dat$kappa = 10
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
dat$lbeta = log(dat$beta)
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[2], lty=1)
#legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[2], lwd=2, lty=rev(1:3), cex=2)
legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=1$,    $\\kappa=0.5$, $w(a_s)\\approx %s$", round(vbGrow(1, 0.5, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[2], lwd=2, lty=rev(1:3), cex=2)

#
#LOGISTIC
#

#Growth
dat$aS = 1 #4 #10
dat$kappa = 0.5 #0.2 #0.1
#Recruitment
dat$alpha = 1.2 
dat$lalpha= log(dat$alpha)
dat$gamma = 1
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.000005
dat$lbeta = log(dat$beta)
#bouncy
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=2, col=cols[1], lty=2, cex.main=3, cex.axis=2, cex.lab=3)
#toned down
dat$aS = 2
dat$kappa = 0.1
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.000005
dat$lbeta = log(dat$beta)
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[1], lty=3)
#prod
dat$aS = 0.1
dat$kappa = 10
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.0000035
dat$lbeta = log(dat$beta)
dat$iterate()
par(mar=c(5, 5, 4, 2)+0.1)
dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[1], lty=1)
legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=1$,    $\\kappa=0.5$, $w(a_s)\\approx %s$", round(vbGrow(1, 0.5, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[1], lwd=2, lty=rev(1:3), cex=2)

#
dev.off()

#
#RP
#

#
xis = seq(0.05, 4, 0.1)

#
png('rpTriptic.png')#, width=480*3)
#layout(t(1:3))

#BH

#
gamma = -1

##
#a = 10
#k = 0.1
#
#lines(xis, getZetaG(xis, M, 0.2, WW, 4, -1, gamma), lty=1, lwd=2)#, col=cols[3])
a = 1 #4
k = 0.5 #0.2
#	        x,   M, k, W, aS, a0, G
line = getZetaG(xis, M, k, WW, a, -1, gamma)
plot(xis, line, 'l', lty=2, lwd=2, ylim=c(0.1, 0.6), col=cols[3],
	xlab = TeX("$F_{MSY}/M$"),
        ylab = TeX('$B_{MSY}/B_0$'),
        main="Space of Reference Points"
)
#
a = 2
k = 0.1
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=3, lwd=2, col=cols[3])
#
a = 0.1
k = 10
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=1, lwd=2, col=cols[3])

#Ricker

#NOTE: there is something going on as gamma passes through 0 and the effect that growth has on RPs.
#inferring with a ricker model has minimal (is any; numerical only? i dont think so.?.) RP effect.
#at some positive gamma value slower growth results in more B^*/B0 not less.  
gamma = 0.0001 

##
#a = 10
#k = 0.1
#
#lines(xis, getZetaG(xis, M, 0.2, WW, 4, -1, gamma), lty=1, lwd=2)#, col=cols[3])
a = 1 #4
k = 0.2 #0.2
#	        x,   M, k, W, aS, a0, G
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=2, lwd=2, col=cols[2])
#
a = 2
k = 0.1
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=3, lwd=2, col=cols[2])
#
a = 0.1
k = 10
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=1, lwd=2, col=cols[2])

#Logistic

#
gamma = 1

##
#a = 10
#k = 0.1
#
#lines(xis, getZetaG(xis, M, 0.2, WW, 4, -1, gamma), lty=1, lwd=2)#, col=cols[3])
a = 1 #4
k = 0.5 #0.2
#	        x,   M, k, W, aS, a0, G
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=2, lwd=2, col=cols[1])
#
a = 2
k = 0.1
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=3, lwd=2, col=cols[1])
#
a = 0.1
k = 10
line = getZetaG(xis, M, k, WW, a, -1, gamma)
lines(xis, line, lty=1, lwd=2, col=cols[1])                         #TeX("$a_s=10$  $\\kappa=0.1$")
legend('bottomleft', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=1$,    $\\kappa=0.5$, $w(a_s)\\approx %s$", round(vbGrow(1, 0.5, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), 
	col=c(rep("black", 3)), 
	#fill=c(cols[1:3]),#, rep(NA, 3)),
	lwd=2, #c(rep(1,3), rep(2,3)), 
	lty=c(rev(1:3)), 
	cex=1
)
legend('top', legend=c("Logistic", "Ricker", "Beverton-Holt"), 
	col=c(cols[1:3]), 
	#fill=c(cols[1:3]),#, rep(NA, 3)),
	lwd=2, #c(rep(1,3), rep(2,3)), 
	lty=c(1, 1, 1), 
	cex=1,
	horiz=T
)
#egend('bottomleft', legend=c("Logistic", "Ricker", "Beverton-Holt", TeX(sprintf("$a_s=2$    $\\kappa=0.1$ $w(a_s)=%s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$   $\\kappa=0.2$ $w(a_s)=%s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$ $\\kappa=10$ $w(a_s)=%s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), 
#       col=c(cols[1:3], rep("black", 3)), 
#       #fill=c(cols[1:3]),#, rep(NA, 3)),
#       lwd=2, #c(rep(1,3), rep(2,3)), 
#       lty=c(1, 1, 1, rev(1:3)), 
#       cex=1
#

dev.off()


#No Impact of growth

#
gStart = 0.0001
lg = function(g){
	#
	a = 0.1
	k = 10
	#               x,   M, k, W, aS, a0, G
	line = getZetaG(xis, M, k, WW, a, -1, g)
	#
	a = 2
	k = 0.1
	line1 = getZetaG(xis, M, k, WW, a, -1, g)
	#
	a = 4
	k = 0.2
	line2 = getZetaG(xis, M, k, WW, a, -1, g)
	#
	return( sum(2*line-line1-line2) )
}
out = uniroot(lg, c(gStart, 1))

