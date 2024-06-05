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
        FF = 0 #FMsy(dat$M, dat$kappa, ww, mod$WW, alpha, beta, gamma) #0 #dat$FMsy
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
        #N = R/(mod$M+FF)
        #Y = ww*R + mod$kappa*(mod$WW*N-B) - (M+FF)*B
        #Y = ww*R + mod$kappa*(mod$WW*R/(M+FF)-B) - (M+FF)*B
        #Y = ww*R + mod$kappa*mod$WW*N - mod$kappa*B - mod$M*B - FF*B
        #fDat = function(x){ wwDat*SRR(x, dat) -dat$M*x-dat$kappa*x + 
        #       dat$kappa*dat$WW*NCar(x, 0, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma) }
        # (alpha*x*( 1-beta*gamma*x )^(1/gamma))/(M+FF)
        Y = R-B*dat$M*(dat$M+dat$kappa)/dat$kappa/dat$WW/(1+dat$M*ww/dat$kappa/dat$WW)
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

##
#png('growthTriptic.png', width=480*3)
#layout(t(1:3))
#
##
##BH
##
#
##Growth
#dat$aS = 4 #10
#dat$kappa = 0.2 #0.1
##Recruitment
#dat$alpha = 1.2
#dat$lalpha= log(dat$alpha)
##dat$beta  = 0.0002
##dat$lbeta = log(dat$beta)
#dat$gamma = -1
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.00006 #-0.00007
#dat$lbeta = log(dat$beta)
#dat$iterate()
#BB = dat$B[46]
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Beverton-Holt", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=2, col=cols[3], cex.main=3, cex.axis=2, cex.lab=3, lty=2)
##toned down
#dat$aS = 2
#dat$kappa = 0.1
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
#dat$lbeta = log(dat$beta)
#dat$iterate()
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[3], lty=3)
##prod
#dat$aS = 0.1
#dat$kappa = 10
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
#dat$lbeta = log(dat$beta)
#dat$iterate()
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[3], lty=1)
#legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[3], lwd=2, lty=rev(1:3), cex=2)
#
##
##RICKER
##
#
##Growth
#dat$aS = 4 #10
#dat$kappa = 0.2 #0.1
##Recruitment
#dat$alpha = 1.2 
#dat$lalpha= log(dat$alpha)
##dat$beta  = 0.0002 
##dat$lbeta = log(dat$beta)
#dat$gamma = 0.0001
#dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)#-0.000015
#dat$lbeta = log(dat$beta)
#dat$iterate()
#BB = dat$B[46]
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Ricker", ylim=c(0,B0), xlab="Time", ylab="Biomass", col=cols[2], lwd=2, cex.main=3, cex.axis=2, cex.lab=3, lty=2)
##toned down
#dat$aS = 2
#dat$kappa = 0.1
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
#dat$lbeta = log(dat$beta)
#dat$iterate()
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[2], lty=3)
##prod
#dat$aS = 0.1
#dat$kappa = 10
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBetaBmsy(BB, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma) #getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)
#dat$lbeta = log(dat$beta)
#dat$iterate()
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[2], lty=1)
#legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[2], lwd=2, lty=rev(1:3), cex=2)

#
#LOGISTIC
#



#Growth
dat$aS = 10 	#4 #10
dat$kappa = 0.1 #0.1
#Recruitment
dat$alpha = 1.2 
dat$lalpha= log(dat$alpha)
dat$gamma = 1
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.000005
dat$lbeta = log(dat$beta)
#
Fmsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
BMsy = BBar(Fmsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)


#
png("shockBiomass.png")

cols = brewer.pal(4, 'Set1')[c(1,2,1,3)]

#
layout(
matrix(c(
	1,1,3,
	1,1,3,
	2,2,4
), 3, 3, byrow=T)
)

##
#dat$time = 1:40
#dat$iterate()
#plot(dat$time-1, dat$B, 'l', col=cols[1], lwd=3, ylim=c(3000,B0), xlab="Time", ylab="Biomass")
##dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[1]) 
##, cex.main=3, cex.axis=2, cex.lab=3)
#points(40, dat$B[40], pch=19, col=cols[1])

#
par(mar=c(4, 4, 4, 2)+0.1)

#
dat$time = 1:31
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[2], add=T)  
#lines(dat$time-1, dat$B, 'l', col=cols[2], lwd=3)
plot(dat$time-1, dat$B, 'l', col=cols[2], lwd=3, ylim=c(3000,B0), ylab="Biomass", xlab="")
points(30, dat$B[30], pch=19, col=cols[2])
text(which(dat$B==max(dat$B[20:30]))-2, max(dat$B[20:30])+200, col=cols[2], label="V", font=2)
#
#points(30, dat$B[30], col=cols[1])

##
#lines((dat$time-1)[21:31]+0.5, SRR(dat$B[11:21], dat), 'l', col=cols[3], lwd=3)
##
#points(30+0.5, SRR(dat$B[21], dat), col=cols[3], pch=19) 
#points(30+0.5, SRR(dat$B[21], dat), col=cols[2], lwd=2)

#
dat$time = 1:21
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[3], add=T)  
#points(11, dat$B[11], col=cols[3])
lines(dat$time-1, dat$B, 'l', col=cols[3], lwd=3)
points(20.5, dat$B[21], pch=19, col=cols[3])
text(which(dat$B==max(dat$B[10:20]))-1, max(dat$B[10:20])+200, col=cols[3], label="V", font=2)
#
points(20.5, dat$B[21], col=cols[2], lwd=2)

##
#lines((dat$time-1)[11:21]+0.5, c(SRR(dat$B[1:11], dat)), 'l', col=cols[4], lwd=3)
#points(10, dat$B[11], col=cols[4], lwd=2)
#points(20+0.5, SRR(dat$B[11], dat), col=cols[4], pch=19)
#points(20+0.5, SRR(dat$B[11], dat), col=cols[3], lwd=2) 

#
dat$time = 1:11
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[4], add=T)  
lines(dat$time-1, dat$B, 'l', col=cols[4], lwd=3)
points(0, dat$B[1], pch=19)
points(0, dat$B[1], col=cols[4], lwd=2)
points(10, dat$B[10], pch=19, col=cols[4])
#
points(10, dat$B[10], col=cols[3], lwd=2)
#
#abline(h=BMsy)
#
lines(dat$time-0.5, SRR(rep(dat$B[1], 11), dat))
points(10+0.5, SRR(dat$B[1], dat), pch=19)
points(10+0.5, SRR(dat$B[1], dat), col=cols[4], lwd=2)

#
dat$time = 1:31
dat$iterate()
#
par(mar=c(4, 4, 0, 2)+0.1)
plot((dat$time-1)[21:31]+0.5, surplus(dat$B[11:21], dat), 'l', col=cols[3], lwd=3, xlab="Time", ylab="Surplus Recruitment", xlim=c(0, 30), ylim=c(0, 2500))#, ylim=c(2250, 3750))
#
points(30+0.5, surplus(dat$B[21], dat), col=cols[3], pch=19) 
points(30+0.5, surplus(dat$B[21], dat), col=cols[2], lwd=2)

#
dat$time = 1:21
dat$iterate()
#
lines((dat$time-1)[11:21]+0.5, c(surplus(dat$B[1:11], dat)), 'l', col=cols[4], lwd=3)
#points(10, dat$B[11], col=cols[4], lwd=2)
points(20+0.5, surplus(dat$B[11], dat), col=cols[4], pch=19)
points(20+0.5, surplus(dat$B[11], dat), col=cols[3], lwd=2) 

#
dat$time = 1:11
dat$iterate()
#
lines(dat$time-0.5, surplus(rep(dat$B[1], 11), dat))
points(10+0.5, surplus(dat$B[1], dat), pch=19)
points(10+0.5, surplus(dat$B[1], dat), col=cols[4], lwd=2)

##
#lines((dat$time-1)[21:31]+0.5, SRR(dat$B[11:21], dat), 'l', col=cols[3], lwd=3)
##
#points(30+0.5, SRR(dat$B[21], dat), col=cols[3], pch=19) 
#points(30+0.5, SRR(dat$B[21], dat), col=cols[2], lwd=2)


#
par(mar=c(4, 4, 4, 2)+0.1)

#
Fmsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)

#
dat$time = 1:44
dat$iterate()
#
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
BMsy = BBar(Fmsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
#
g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
g = Vectorize(g, 'x')
maxF = uniroot(g, c(Fmsy, exp(dat$lalpha)))$root
FFs = seq(0, maxF, length.out=1000)
#
#plot(FFs*g(FFs), g(FFs), type='l', lwd=2, xlab="Yield", ylab="Biomass", ylim=c(0, 10000))
BB = seq(0, 15000, length.out=1000)
base = 0
plot(surplus(BB, dat), BB, type='l', lwd=2, xlab="Surplus Recruitment", ylab="Biomass", ylim=c(3000, 10000), xlim=c(base, 2500))#max(surplus(BB, dat))) )
#
segments(base+0, dat$B[1], base+0, dat$B[11], col=cols[4], lwd=3)
points(base+0, dat$B[1], pch=19, col='black')
points(base+0, dat$B[1], col=cols[4], lwd=2)
points(base+0, dat$B[11], pch=19, col=cols[4])
segments(base+100, dat$B[11], base+100, max(dat$B[11:22]), col=cols[3], lwd=3)
points(base+100, dat$B[11], col="white", pch=19)
points(base+100, dat$B[11], col=cols[3], lwd=2)
text(base+100, max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
points(base+100, dat$B[22], pch=19, col=cols[3])
segments(base+200, dat$B[22], base+200, max(dat$B[22:33]), col=cols[2], lwd=3)
points(base+200, dat$B[22], col="white", pch=19)
points(base+200, dat$B[22], col=cols[2], lwd=2)
text(base+200, max(dat$B[22:33])+150, label="V", col=cols[2], font=2)
points(base+200, dat$B[33], pch=19, col=cols[2])
#segments(300, dat$B[33], 300, max(dat$B[33:44]), col=cols[1], lwd=3)
#points(300, dat$B[33], col=cols[1])
#points(300, max(dat$B[33:44]), pch="V", col=cols[1])
#points(300, dat$B[44], pch=19, col=cols[1])
#abline(h=BMsy)
#
#BB = seq(dat$B[1], dat$B[11], length.out=1000)
#lines(surplus(BB, dat), BB, col=cols[3], lwd=3)
points(surplus(dat$B[11], dat), dat$B[11], col=cols[4], pch=19)
points(surplus(dat$B[11], dat), dat$B[11], col=cols[3], lwd=2)
#
#BB = seq(dat$B[11], dat$B[22], length.out=1000)
points(surplus(dat$B[22], dat), dat$B[22], col=cols[3], pch=19)
points(surplus(dat$B[22], dat), dat$B[22], col=cols[2], lwd=2)
#points(surplus(max(dat$B[11:22])+150, dat), max(dat$B[11:22])+150, pch="V", col=cols[3], lwd=2)
text(surplus(max(dat$B[11:22])+150, dat), max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
#
points(surplus(dat$B[33], dat), dat$B[33], col=cols[2], pch=19)
#points(surplus(max(dat$B[22:33])+250, dat), max(dat$B[22:33])+250, pch="V", col=cols[2], lwd=2)
text(surplus(max(dat$B[22:33])+250, dat), max(dat$B[22:33])+250, label="V", col=cols[2], font=2)
#
#lines(c(0, surplus(BMsy, dat)), c(0, BMsy), col='red')

#
dev.off()

#
png("shockBiomassYeild.png")

cols = brewer.pal(4, 'Set1')[c(1,2,1,3)]

#
layout(
matrix(c(
	1,1,3,
	1,1,3,
	2,2,4
), 3, 3, byrow=T)
)

##
#dat$time = 1:40
#dat$iterate()
#plot(dat$time-1, dat$B, 'l', col=cols[1], lwd=3, ylim=c(3000,B0), xlab="Time", ylab="Biomass")
##dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[1]) 
##, cex.main=3, cex.axis=2, cex.lab=3)
#points(40, dat$B[40], pch=19, col=cols[1])

#
par(mar=c(4, 4, 4, 2)+0.1)

#
dat$time = 1:31
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[2], add=T)  
#lines(dat$time-1, dat$B, 'l', col=cols[2], lwd=3)
plot(dat$time-1, dat$B, 'l', col=cols[2], lwd=3, ylim=c(3000,B0), ylab="Biomass", xlab="")
title("                                                             Oscillation Mechanism")#line = -2)
points(30, dat$B[30], pch=19, col=cols[2])
text(which(dat$B==max(dat$B[20:30]))-2, max(dat$B[20:30])+200, col=cols[2], label="V", font=2)
#
#points(30, dat$B[30], col=cols[1])

##
#lines((dat$time-1)[21:31]+0.5, SRR(dat$B[11:21], dat), 'l', col=cols[3], lwd=3)
##
#points(30+0.5, SRR(dat$B[21], dat), col=cols[3], pch=19) 
#points(30+0.5, SRR(dat$B[21], dat), col=cols[2], lwd=2)

#
dat$time = 1:21
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[3], add=T)  
#points(11, dat$B[11], col=cols[3])
lines(dat$time-1, dat$B, 'l', col=cols[3], lwd=3)
points(20.5, dat$B[21], pch=19, col=cols[3])
text(which(dat$B==max(dat$B[10:20]))-1, max(dat$B[10:20])+200, col=cols[3], label="V", font=2)
#
points(20.5, dat$B[21], col=cols[2], lwd=2)

##
#lines((dat$time-1)[11:21]+0.5, c(SRR(dat$B[1:11], dat)), 'l', col=cols[4], lwd=3)
#points(10, dat$B[11], col=cols[4], lwd=2)
#points(20+0.5, SRR(dat$B[11], dat), col=cols[4], pch=19)
#points(20+0.5, SRR(dat$B[11], dat), col=cols[3], lwd=2) 

#
dat$time = 1:11
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[4], add=T)  
lines(dat$time-1, dat$B, 'l', col=cols[4], lwd=3)
points(0, dat$B[1], pch=19)
points(0, dat$B[1], col=cols[4], lwd=2)
points(10, dat$B[10], pch=19, col=cols[4])
#
points(10, dat$B[10], col=cols[3], lwd=2)
#
#abline(h=BMsy)
#
lines(dat$time-0.5, SRR(rep(dat$B[1], 11), dat))
points(10+0.5, SRR(dat$B[1], dat), pch=19)
points(10+0.5, SRR(dat$B[1], dat), col=cols[4], lwd=2)


#
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
Fmsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
BMsy = BBar(Fmsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
#
r = (Fmsy*BMsy)/surplus(BMsy, dat)

#
dat$time = 1:31
dat$iterate()
#
par(mar=c(4, 4, 0, 2)+0.1)
plot((dat$time-1)[21:31]+0.5, surplus(dat$B[11:21], dat)*r, 'l', col=cols[3], lwd=3, xlab="Time", ylab="Yield", xlim=c(0, 30), ylim=c(0, 1800)) #Fmsy*BMsy))#, ylim=c(2250, 3750))
#
points(30+0.5, surplus(dat$B[21], dat)*r, col=cols[3], pch=19) 
points(30+0.5, surplus(dat$B[21], dat)*r, col=cols[2], lwd=2)

#
dat$time = 1:21
dat$iterate()
#
lines((dat$time-1)[11:21]+0.5, c(surplus(dat$B[1:11], dat)*r), 'l', col=cols[4], lwd=3)
#points(10, dat$B[11], col=cols[4], lwd=2)
points(20+0.5, surplus(dat$B[11], dat)*r, col=cols[4], pch=19)
points(20+0.5, surplus(dat$B[11], dat)*r, col=cols[3], lwd=2) 

#
dat$time = 1:11
dat$iterate()
#
lines(dat$time-0.5, surplus(rep(dat$B[1], 11), dat)*r)
points(10+0.5, surplus(dat$B[1], dat)*r, pch=19)
points(10+0.5, surplus(dat$B[1], dat)*r, col=cols[4], lwd=2)

##
#lines((dat$time-1)[21:31]+0.5, SRR(dat$B[11:21], dat), 'l', col=cols[3], lwd=3)
##
#points(30+0.5, SRR(dat$B[21], dat), col=cols[3], pch=19) 
#points(30+0.5, SRR(dat$B[21], dat), col=cols[2], lwd=2)


#
par(mar=c(4, 4, 4, 2)+0.1)

#
dat$time = 1:44
dat$iterate()

#
g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
g = Vectorize(g, 'x')
maxF = uniroot(g, c(Fmsy, exp(dat$lalpha)))$root
FFs = seq(0, maxF, length.out=1000)
#
#plot(FFs*g(FFs), g(FFs), type='l', lwd=2, xlab="Yield", ylab="Biomass", ylim=c(0, 10000))
BB = seq(0, 15000, length.out=1000)
base = 0
plot(FFs*g(FFs), g(FFs), type='l', lwd=2, xlab="Yield", ylab="Biomass", ylim=c(3000, 10000), xlim=c(base, 1700)) #Fmsy*BMsy))#max(surplus(BB, dat))) )
#
segments(base+0, dat$B[1], base+0, dat$B[11], col=cols[4], lwd=3)
points(base+0, dat$B[1], pch=19, col='black')
points(base+0, dat$B[1], col=cols[4], lwd=2)
points(base+0, dat$B[11], pch=19, col=cols[4])
segments(base+100, dat$B[11], base+100, max(dat$B[11:22]), col=cols[3], lwd=3)
points(base+100, dat$B[11], col="white", pch=19)
points(base+100, dat$B[11], col=cols[3], lwd=2)
text(base+100, max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
points(base+100, dat$B[22], pch=19, col=cols[3])
segments(base+200, dat$B[22], base+200, max(dat$B[22:33]), col=cols[2], lwd=3)
points(base+200, dat$B[22], col="white", pch=19)
points(base+200, dat$B[22], col=cols[2], lwd=2)
text(base+200, max(dat$B[22:33])+150, label="V", col=cols[2], font=2)
points(base+200, dat$B[33], pch=19, col=cols[2])
#segments(300, dat$B[33], 300, max(dat$B[33:44]), col=cols[1], lwd=3)
#points(300, dat$B[33], col=cols[1])
#points(300, max(dat$B[33:44]), pch="V", col=cols[1])
#points(300, dat$B[44], pch=19, col=cols[1])
#abline(h=BMsy)
#
#BB = seq(dat$B[1], dat$B[11], length.out=1000)
#lines(surplus(BB, dat), BB, col=cols[3], lwd=3)
points(surplus(dat$B[11], dat)*r, dat$B[11], col=cols[4], pch=19)
points(surplus(dat$B[11], dat)*r, dat$B[11], col=cols[3], lwd=2)
#
#BB = seq(dat$B[11], dat$B[22], length.out=1000)
points(surplus(dat$B[22], dat)*r, dat$B[22], col=cols[3], pch=19)
points(surplus(dat$B[22], dat)*r, dat$B[22], col=cols[2], lwd=2)
#points(surplus(max(dat$B[11:22])+150, dat), max(dat$B[11:22])+150, pch="V", col=cols[3], lwd=2)
text(surplus(max(dat$B[11:22])+150, dat)*r, max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
#
points(surplus(dat$B[33], dat)*r, dat$B[33], col=cols[2], pch=19)
#points(surplus(max(dat$B[22:33])+250, dat), max(dat$B[22:33])+250, pch="V", col=cols[2], lwd=2)
text(surplus(max(dat$B[22:33])+250, dat)*r, max(dat$B[22:33])+250, label="V", col=cols[2], font=2)
#
#lines(c(0, surplus(BMsy, dat)), c(0, BMsy), col='red')

#
dev.off()


#
png("shockBiomassSurplus.png")

cols = brewer.pal(4, 'Set1')[c(1,2,1,3)]

#
layout(
matrix(c(
	1,1,3,
	1,1,3,
	2,2,4
), 3, 3, byrow=T)
)

##
#dat$time = 1:40
#dat$iterate()
#plot(dat$time-1, dat$B, 'l', col=cols[1], lwd=3, ylim=c(3000,B0), xlab="Time", ylab="Biomass")
##dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[1]) 
##, cex.main=3, cex.axis=2, cex.lab=3)
#points(40, dat$B[40], pch=19, col=cols[1])

#
par(mar=c(4, 4, 4, 2)+0.1)

#
dat$time = 1:31
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[2], add=T)  
#lines(dat$time-1, dat$B, 'l', col=cols[2], lwd=3)
plot(dat$time-1, dat$B, 'l', col=cols[2], lwd=3, ylim=c(3000,B0), ylab="Biomass", xlab="")
title("                                                             Oscillation Mechanism")#line = -2)
points(30, dat$B[30], pch=19, col=cols[2])
text(which(dat$B==max(dat$B[20:30]))-2, max(dat$B[20:30])+200, col=cols[2], label="V", font=2)
#
#points(30, dat$B[30], col=cols[1])

##
#lines((dat$time-1)[21:31]+0.5, SRR(dat$B[11:21], dat), 'l', col=cols[3], lwd=3)
##
#points(30+0.5, SRR(dat$B[21], dat), col=cols[3], pch=19) 
#points(30+0.5, SRR(dat$B[21], dat), col=cols[2], lwd=2)

#
dat$time = 1:21
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[3], add=T)  
#points(11, dat$B[11], col=cols[3])
lines(dat$time-1, dat$B, 'l', col=cols[3], lwd=3)
points(20.5, dat$B[21], pch=19, col=cols[3])
text(which(dat$B==max(dat$B[10:20]))-1, max(dat$B[10:20])+200, col=cols[3], label="V", font=2)
#
points(20.5, dat$B[21], col=cols[2], lwd=2)

##
#lines((dat$time-1)[11:21]+0.5, c(SRR(dat$B[1:11], dat)), 'l', col=cols[4], lwd=3)
#points(10, dat$B[11], col=cols[4], lwd=2)
#points(20+0.5, SRR(dat$B[11], dat), col=cols[4], pch=19)
#points(20+0.5, SRR(dat$B[11], dat), col=cols[3], lwd=2) 

#
dat$time = 1:11
dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[4], add=T)  
lines(dat$time-1, dat$B, 'l', col=cols[4], lwd=3)
points(0, dat$B[1], pch=19)
points(0, dat$B[1], col=cols[4], lwd=2)
points(10, dat$B[10], pch=19, col=cols[4])
#
points(10, dat$B[10], col=cols[3], lwd=2)
#
#abline(h=BMsy)
#
lines(dat$time-0.5, SRR(rep(dat$B[1], 11), dat))
points(10+0.5, SRR(dat$B[1], dat), pch=19)
points(10+0.5, SRR(dat$B[1], dat), col=cols[4], lwd=2)


#
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
Fmsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
BMsy = BBar(Fmsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
#
r = (Fmsy*BMsy)/surplus(BMsy, dat)

#
dat$time = 1:31
dat$iterate()
#
par(mar=c(4, 4, 0, 2)+0.1)
plot((dat$time-1)[21:31]+0.5, surplus(dat$B[11:21], dat)*r, 'l', col=cols[3], lwd=3, xlab="Time", ylab="Surplus Biomass", xlim=c(0, 30), ylim=c(0, 1800)) #Fmsy*BMsy))#, ylim=c(2250, 3750))
#
points(30+0.5, surplus(dat$B[21], dat)*r, col=cols[3], pch=19) 
points(30+0.5, surplus(dat$B[21], dat)*r, col=cols[2], lwd=2)

#
dat$time = 1:21
dat$iterate()
#
lines((dat$time-1)[11:21]+0.5, c(surplus(dat$B[1:11], dat)*r), 'l', col=cols[4], lwd=3)
#points(10, dat$B[11], col=cols[4], lwd=2)
points(20+0.5, surplus(dat$B[11], dat)*r, col=cols[4], pch=19)
points(20+0.5, surplus(dat$B[11], dat)*r, col=cols[3], lwd=2) 

#
dat$time = 1:11
dat$iterate()
#
lines(dat$time-0.5, surplus(rep(dat$B[1], 11), dat)*r)
points(10+0.5, surplus(dat$B[1], dat)*r, pch=19)
points(10+0.5, surplus(dat$B[1], dat)*r, col=cols[4], lwd=2)

##
#lines((dat$time-1)[21:31]+0.5, SRR(dat$B[11:21], dat), 'l', col=cols[3], lwd=3)
##
#points(30+0.5, SRR(dat$B[21], dat), col=cols[3], pch=19) 
#points(30+0.5, SRR(dat$B[21], dat), col=cols[2], lwd=2)


#
par(mar=c(4, 4, 4, 2)+0.1)

#
dat$time = 1:44
dat$iterate()

#
g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
g = Vectorize(g, 'x')
maxF = uniroot(g, c(Fmsy, exp(dat$lalpha)))$root
FFs = seq(0, maxF, length.out=1000)
#
#plot(FFs*g(FFs), g(FFs), type='l', lwd=2, xlab="Yield", ylab="Biomass", ylim=c(0, 10000))
BB = seq(0, 15000, length.out=1000)
base = 0
plot(FFs*g(FFs), g(FFs), type='l', lwd=2, xlab="Surplus Biomass", ylab="Biomass", ylim=c(3000, 10000), xlim=c(base, 1700)) #Fmsy*BMsy))#max(surplus(BB, dat))) )
#
segments(base+0, dat$B[1], base+0, dat$B[11], col=cols[4], lwd=3)
points(base+0, dat$B[1], pch=19, col='black')
points(base+0, dat$B[1], col=cols[4], lwd=2)
points(base+0, dat$B[11], pch=19, col=cols[4])
segments(base+100, dat$B[11], base+100, max(dat$B[11:22]), col=cols[3], lwd=3)
points(base+100, dat$B[11], col="white", pch=19)
points(base+100, dat$B[11], col=cols[3], lwd=2)
text(base+100, max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
points(base+100, dat$B[22], pch=19, col=cols[3])
segments(base+200, dat$B[22], base+200, max(dat$B[22:33]), col=cols[2], lwd=3)
points(base+200, dat$B[22], col="white", pch=19)
points(base+200, dat$B[22], col=cols[2], lwd=2)
text(base+200, max(dat$B[22:33])+150, label="V", col=cols[2], font=2)
points(base+200, dat$B[33], pch=19, col=cols[2])
#segments(300, dat$B[33], 300, max(dat$B[33:44]), col=cols[1], lwd=3)
#points(300, dat$B[33], col=cols[1])
#points(300, max(dat$B[33:44]), pch="V", col=cols[1])
#points(300, dat$B[44], pch=19, col=cols[1])
#abline(h=BMsy)
#
#BB = seq(dat$B[1], dat$B[11], length.out=1000)
#lines(surplus(BB, dat), BB, col=cols[3], lwd=3)
points(surplus(dat$B[11], dat)*r, dat$B[11], col=cols[4], pch=19)
points(surplus(dat$B[11], dat)*r, dat$B[11], col=cols[3], lwd=2)
#
#BB = seq(dat$B[11], dat$B[22], length.out=1000)
points(surplus(dat$B[22], dat)*r, dat$B[22], col=cols[3], pch=19)
points(surplus(dat$B[22], dat)*r, dat$B[22], col=cols[2], lwd=2)
#points(surplus(max(dat$B[11:22])+150, dat), max(dat$B[11:22])+150, pch="V", col=cols[3], lwd=2)
text(surplus(max(dat$B[11:22])+150, dat)*r, max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
#
points(surplus(dat$B[33], dat)*r, dat$B[33], col=cols[2], pch=19)
#points(surplus(max(dat$B[22:33])+250, dat), max(dat$B[22:33])+250, pch="V", col=cols[2], lwd=2)
text(surplus(max(dat$B[22:33])+250, dat)*r, max(dat$B[22:33])+250, label="V", col=cols[2], font=2)
#
#lines(c(0, surplus(BMsy, dat)), c(0, BMsy), col='red')

#
dev.off()







#toned down
#dat$aS = 2
#dat$kappa = 0.1
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.000005
#dat$lbeta = log(dat$beta)
#dat$iterate()
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[1], lty=3)
##prod
#dat$aS = 0.1
#dat$kappa = 10
#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
#dat$beta = getBeta(B0, M, dat$kappa, ww, dat$WW, dat$alpha, dat$gamma)+0.0000035
#dat$lbeta = log(dat$beta)
#dat$iterate()
#par(mar=c(5, 5, 4, 2)+0.1)
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass", add=T, lwd=2, col=cols[1], lty=1)
#legend('bottom', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), col=cols[1], lwd=2, lty=rev(1:3), cex=2)

#
png("shockYield.png")
#par(mar=c(5, 5, 4, 2)+0.1)

#
dat$time = 1:44
dat$iterate()
#
g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
g = Vectorize(g, 'x')
maxF = uniroot(g, c(Fmsy, exp(dat$lalpha)))$root
FFs = seq(0, maxF, length.out=1000)
#
plot(g(FFs), FFs*g(FFs), type='l', lwd=2, xlab="Biomass", ylab="Equilibrium Surplus Biomass", main="Yield Curve")
#
segments(dat$B[1], 0, dat$B[11], 0, col=cols[4], lwd=3)
points(dat$B[1], 0, col=cols[4])
points(dat$B[11], 0, pch=19, col=cols[4])
segments(dat$B[11], 100, max(dat$B[11:22]), 100, col=cols[3], lwd=3)
points(dat$B[11], 100, col=cols[3])
points(max(dat$B[11:22]), 100, pch="<", col=cols[3])
points(dat$B[22], 100, pch=19, col=cols[3])
segments(dat$B[22], 200, max(dat$B[22:33]), 200, col=cols[2], lwd=3)
points(dat$B[22], 200, col=cols[2])
points(max(dat$B[22:33]), 200, pch="<", col=cols[2])
points(dat$B[33], 200, pch=19, col=cols[2])
#segments(dat$B[33], 300, max(dat$B[33:44]), 300, col=cols[1], lwd=3)
#points(dat$B[33], 300, col=cols[1])
#points(max(dat$B[33:44]), 300, pch="<", col=cols[1])
#points(dat$B[44], 300, pch=19, col=cols[1])

#
dev.off()

#
png("shockYieldInv.png")
#par(mar=c(5, 5, 4, 2)+0.1)

#
Fmsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)

#
dat$time = 1:44
dat$iterate()
#
ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
BMsy = BBar(Fmsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
#
g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
g = Vectorize(g, 'x')
maxF = uniroot(g, c(Fmsy, exp(dat$lalpha)))$root
FFs = seq(0, maxF, length.out=1000)
#
#plot(FFs*g(FFs), g(FFs), type='l', lwd=2, xlab="Yield", ylab="Biomass", ylim=c(0, 10000))
BB = seq(0, max(dat$B), length.out=1000)
plot(surplus(BB, dat), BB, type='l', lwd=2, xlab=TeX("\\frac{dB}{dt}"), ylab="Biomass", ylim=c(0, 10000))
#
segments(0, dat$B[1], 0, dat$B[11], col=cols[4], lwd=3)
points(0, dat$B[1], col=cols[4], lwd=2)
points(0, dat$B[11], pch=19, col=cols[4])
segments(100, dat$B[11], 100, max(dat$B[11:22]), col=cols[3], lwd=3)
points(100, dat$B[11], col=cols[3], lwd=2)
text(100, max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
points(100, dat$B[22], pch=19, col=cols[3])
segments(200, dat$B[22], 200, max(dat$B[22:33]), col=cols[2], lwd=3)
points(200, dat$B[22], col=cols[2], lwd=2)
text(200, max(dat$B[22:33])+150, label="V", col=cols[2], font=2)
points(200, dat$B[33], pch=19, col=cols[2])
#segments(300, dat$B[33], 300, max(dat$B[33:44]), col=cols[1], lwd=3)
#points(300, dat$B[33], col=cols[1])
#points(300, max(dat$B[33:44]), pch="V", col=cols[1])
#points(300, dat$B[44], pch=19, col=cols[1])
#abline(h=BMsy)
#
#BB = seq(dat$B[1], dat$B[11], length.out=1000)
#lines(surplus(BB, dat), BB, col=cols[3], lwd=3)
points(surplus(dat$B[11], dat), dat$B[11], col=cols[4], pch=19)
points(surplus(dat$B[11], dat), dat$B[11], col=cols[3], lwd=2)
#
#BB = seq(dat$B[11], dat$B[22], length.out=1000)
points(surplus(dat$B[22], dat), dat$B[22], col=cols[3], pch=19)
points(surplus(dat$B[22], dat), dat$B[22], col=cols[2], lwd=2)
#points(surplus(max(dat$B[11:22])+150, dat), max(dat$B[11:22])+150, pch="V", col=cols[3], lwd=2)
text(surplus(max(dat$B[11:22])+150, dat), max(dat$B[11:22])+150, label="V", col=cols[3], font=2)
#
points(surplus(dat$B[33], dat), dat$B[33], col=cols[2], pch=19)
#points(surplus(max(dat$B[22:33])+250, dat), max(dat$B[22:33])+250, pch="V", col=cols[2], lwd=2)
text(surplus(max(dat$B[22:33])+250, dat), max(dat$B[22:33])+250, label="V", col=cols[2], font=2)
#
lines(c(0, surplus(BMsy, dat)), c(0, BMsy), col='red')

#
dev.off()

##
#m1F = 0
#m2F = uniroot(function(x){g(x)-dat$B[11]}, c(Fmsy, exp(dat$lalpha)))$root
#FFs = seq(m1F, m2F, length.out=1000)
#lines(g(FFs), FFs*g(FFs), col=cols[4], lwd=3)
##points(g(FFs)[1], (FFs*g(FFs))[1], pch="|", col=cols[4])
#points(g(FFs)[1000], (FFs*g(FFs))[1000], pch=19, col=cols[4])
##
##F = uniroot(function(x){g(x)-dat$B[12]}, c(Fmsy, exp(dat$lalpha)))$root
#m1F = m2F
#m2F = uniroot(function(x){g(x)-max(dat$B[12:22])}, c(0, exp(dat$lalpha)))$root
##m3F = uniroot(function(x){g(x)-dat$B[22]}, c(0, exp(dat$lalpha)))$root
#FFs = seq(m1F, m2F, length.out=1000)
#lines(g(FFs), FFs*g(FFs), col=cols[3], lwd=3)
#m3F = uniroot(function(x){g(x)-dat$B[22]}, c(0, exp(dat$lalpha)))$root
#FFs = seq(m2F, m3F, length.out=1000)
#points(g(FFs)[1000], (FFs*g(FFs))[1000], pch=19, col=cols[3])
##points(g(FFs)[1000], (FFs*g(FFs))[1000], pch="<", col=cols[3])
##
#m1F = m3F
#m2F = uniroot(function(x){g(x)-max(dat$B[23:33])}, c(0, exp(dat$lalpha)))$root
##m3F = uniroot(function(x){g(x)-dat$B[22]}, c(0, exp(dat$lalpha)))$root
#FFs = seq(m1F, m2F, length.out=1000)
#lines(g(FFs), FFs*g(FFs), col=cols[2], lwd=3)
#m3F = uniroot(function(x){g(x)-dat$B[33]}, c(0, exp(dat$lalpha)))$root
#FFs = seq(m2F, m3F, length.out=1000)
#points(g(FFs)[1000], (FFs*g(FFs))[1000], pch=19, col=cols[2])
##
#m1F = m3F
#m2F = uniroot(function(x){g(x)-max(dat$B[34:44])}, c(0, exp(dat$lalpha)))$root
##m3F = uniroot(function(x){g(x)-dat$B[22]}, c(0, exp(dat$lalpha)))$root
#FFs = seq(m1F, m2F, length.out=1000)
#lines(g(FFs), FFs*g(FFs), col=cols[1], lwd=3)
#m3F = uniroot(function(x){g(x)-dat$B[44]}, c(0, exp(dat$lalpha)))$root
#FFs = seq(m2F, m3F, length.out=1000)
#points(g(FFs)[1000], (FFs*g(FFs))[1000], pch=19, col=cols[1])
##
#abline(v=g(Fmsy))



#BBar = (FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#lines(dat$B[1:11],  col=cols[4], lwd=3)
#segments(dat$B[1], 0, dat$B[11], 0, col=cols[4], lwd=3)
#segments(dat$B[12], 0, max(dat$B[12:22]), 0, col=cols[3], lwd=3)
#segments(dat$B[23], 0, dat$B[33], 0, col=cols[2], lwd=3)
#segments(dat$B[34], 0, dat$B[44], 0, col=cols[1], lwd=3)
#rug(dat$B[34:44], col=cols[1], lwd=3)
#rug(dat$B[23:33], col=cols[2], lwd=3)
#rug(dat$B[12:22], col=cols[3], lwd=3)
#rug(dat$B[1:11], col=cols[4], lwd=3)

##
#dat$time = 1:33
#dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[2], add=T)
#
##
#dat$time = 1:22
#dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[3], add=T)
#
##
#dat$time = 1:11
#dat$iterate()
#dat$plotQuan( function(B){B}, main="Logistic", ylim=c(0,B0), xlab="Time", ylab="Biomass", lwd=3, col=cols[4], add=T)
#

#
#dev.off()

##
##RP
##
#
##
#xis = seq(0.05, 4, 0.1)
#
##
#png('rpTriptic.png')#, width=480*3)
##layout(t(1:3))
#
##BH
#
##
#gamma = -1
#
###
##a = 10
##k = 0.1
##
##lines(xis, getZetaG(xis, M, 0.2, WW, 4, -1, gamma), lty=1, lwd=2)#, col=cols[3])
#a = 4
#k = 0.2
##	        x,   M, k, W, aS, a0, G
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#plot(xis, line, 'l', lty=2, lwd=2, ylim=c(0.1, 0.6), col=cols[3],
#	xlab = TeX("$F_{MSY}/M$"),
#        ylab = TeX('$B_{MSY}/B_0$'),
#        main="Space of Reference Points"
#)
##
#a = 2
#k = 0.1
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=3, lwd=2, col=cols[3])
##
#a = 0.1
#k = 10
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=1, lwd=2, col=cols[3])
#
##Ricker
#
##
#gamma = 0.0001
#
###
##a = 10
##k = 0.1
##
##lines(xis, getZetaG(xis, M, 0.2, WW, 4, -1, gamma), lty=1, lwd=2)#, col=cols[3])
#a = 4
#k = 0.2
##	        x,   M, k, W, aS, a0, G
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=2, lwd=2, col=cols[2])
##
#a = 2
#k = 0.1
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=3, lwd=2, col=cols[2])
##
#a = 0.1
#k = 10
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=1, lwd=2, col=cols[2])
#
##Logistic
#
##
#gamma = 1
#
###
##a = 10
##k = 0.1
##
##lines(xis, getZetaG(xis, M, 0.2, WW, 4, -1, gamma), lty=1, lwd=2)#, col=cols[3])
#a = 4
#k = 0.2
##	        x,   M, k, W, aS, a0, G
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=2, lwd=2, col=cols[1])
##
#a = 2
#k = 0.1
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=3, lwd=2, col=cols[1])
##
#a = 0.1
#k = 10
#line = getZetaG(xis, M, k, WW, a, -1, gamma)
#lines(xis, line, lty=1, lwd=2, col=cols[1])                         #TeX("$a_s=10$  $\\kappa=0.1$")
#legend('bottomleft', legend=c(TeX(sprintf("$a_s=2$,    $\\kappa=0.1$, $w(a_s)\\approx %s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$,    $\\kappa=0.2$, $w(a_s)\\approx %s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$, $\\kappa=10$,  $w(a_s)\\approx %s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), 
#	col=c(rep("black", 3)), 
#	#fill=c(cols[1:3]),#, rep(NA, 3)),
#	lwd=2, #c(rep(1,3), rep(2,3)), 
#	lty=c(rev(1:3)), 
#	cex=1
#)
#legend('top', legend=c("Logistic", "Ricker", "Beverton-Holt"), 
#	col=c(cols[1:3]), 
#	#fill=c(cols[1:3]),#, rep(NA, 3)),
#	lwd=2, #c(rep(1,3), rep(2,3)), 
#	lty=c(1, 1, 1), 
#	cex=1,
#	horiz=T
#)
##egend('bottomleft', legend=c("Logistic", "Ricker", "Beverton-Holt", TeX(sprintf("$a_s=2$    $\\kappa=0.1$ $w(a_s)=%s$", round(vbGrow(2, 0.1, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=4$   $\\kappa=0.2$ $w(a_s)=%s$", round(vbGrow(4, 0.2, dat$WW, dat$a0), 2))), TeX(sprintf("$a_s=0.1$ $\\kappa=10$ $w(a_s)=%s$", round(vbGrow(0.1, 10, dat$WW, dat$a0),2)))), 
##       col=c(cols[1:3], rep("black", 3)), 
##       #fill=c(cols[1:3]),#, rep(NA, 3)),
##       lwd=2, #c(rep(1,3), rep(2,3)), 
##       lty=c(1, 1, 1, rev(1:3)), 
##       cex=1
##
#
#dev.off()
#
#
##No Impact of growth
#
##
#gStart = 0.0001
#lg = function(g){
#	#
#	a = 0.1
#	k = 10
#	#               x,   M, k, W, aS, a0, G
#	line = getZetaG(xis, M, k, WW, a, -1, g)
#	#
#	a = 2
#	k = 0.1
#	line1 = getZetaG(xis, M, k, WW, a, -1, g)
#	#
#	a = 4
#	k = 0.2
#	line2 = getZetaG(xis, M, k, WW, a, -1, g)
#	#
#	return( sum(2*line-line1-line2) )
#}
#out = uniroot(lg, c(gStart, 1))

