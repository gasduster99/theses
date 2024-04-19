rm(list=ls())

#
library(VGAM)
library(boot)
library(pracma)
library(mvtnorm)
library(graphics)
library(parallel)
library(latex2exp)
library(rootSolve)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)
#
#source('gpFunk.r')
source('gpFunkGaussLxLy.r')

#
#FUNCTIONS
#

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
SRR = function(B, alpha, beta, gamma){
        alpha*B*(1-beta*gamma*B)^(1/gamma)
}

#
PBar = function(alpha, beta, gamma, ff, M){ 1/(gamma*beta)*(1-((M+ff)/alpha)^gamma) } #((alpha/(M+

#
getBeta = function(alpha, gamma, M, B0){
        (1-(M/alpha)^gamma)/B0/gamma
}
b = Vectorize(getBeta, c("alpha", "gamma"))

#
getAlpha = function(gamma, ff, M){
        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}
a = Vectorize(getAlpha, "gamma")

#
getZeta = function(gamma, ff, M){
        #fgfm = ff*gamma/(ff+M)
        #fgfm/(1 + fgfm - (M/ff+M)^gamma)

        #
        (1-((M+ff)/getAlpha(gamma, ff, M))^gamma) / (1-(M/getAlpha(gamma, ff, M))^gamma)
}
z = Vectorize(getZeta, "gamma", "ff")

#
FMsy = function(alpha, gamma, M){
        #
	if(gamma==-1){ return(sqrt(alpha*M)-M) }
	#
        FUpper = alpha-M
        root = uniroot.all(function(ff){ a(gamma, ff, M)-alpha }, c(0, FUpper))

        #root = uniroot.all(function(ff){ 1 - ff/gamma/(M+ff) - (alpha/(M+ff))^(-1/gamma) }, c(0, 
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        ##
        return(root)
}

#
FSPRX = function(M, x){ 
	M*((1/x)-1)
}

#
getAlphaProxy = function(gamma, M, x, y){
	M*( (1/(x^gamma)-y)/(1-y) )^(1/gamma)
}

#
getAlphaMsy = function(gamma, ff, M){
        (ff+M)*(1+ff*gamma/(ff+M))^(1/gamma)
}

#
getGammaProMsy = function(x, y){
	pow = y/(x-1)/(y-1)
	pow - lambertWn(pow*x^pow*log(x))/log(x)
}

#
#HEAD
#

#
names = c("RF", "GF", "FF")
xs = c(0.50, 0.45, 0.30)
ys = c(0.40, 0.40, 0.25)
#RF, GF, Flatfish
who = 1 #3
who = 1 #3
name = names[who]
x = xs[who]
y = ys[who]
#
M = 0.2
#0.2, 0.2444..., 0.4666...
FSPRx = M*((1/x)-1) 

#
gs = seq(-1.3, 0.1, length.out=1000)#seq(-2, 2, length.out=100)
aMsy = getAlpha(gs, FSPRx, M)

#alphas and gammas: to hit proxy
#alphas and gammas: to hit MSY
#point (alpha, gamma): to hit both proxy and MSY

#
aProxy = getAlphaProxy(gs, M, x, y)
#
gW = getGammaProMsy(x, y)

#
png(sprintf("ag%s_x%s_y%s.png", name, x, y))
plot(gs, aMsy, 'l', 
	lwd=3, 
	ylim=c(min(aProxy), 6*M),
	xlab=TeX("$\\gamma$"),
	ylab=TeX("$\\alpha$"),
	main=TeX("Schnute $\\alpha$-$\\gamma$ Relationships")
)
lines(gs, aProxy, lwd=3, lty=2)
points(gW, getAlpha(gW, FSPRx, M), col='red', pch=19)
#abline(v=-1)
points(-1, getAlpha(-1, FSPRx, M), col='blue', pch=19)
points(-1, getAlphaProxy(-1, M, x, y), col='blue', pch=19)
segments(-1, getAlpha(-1, FSPRx, M), -1, getAlphaProxy(-1, M, x, y), col='blue')
legend("topright", 
	legend=c("MSY", "Proxy", "Proxy=MSY", "BH"), 
	lty=c(1, 2, NA, NA), 
	lwd=3, 
	pch=c(NA, NA, 19, 19),
	col=c('black', 'black', 'red', 'blue')
)
#points(gW, getAlpha(gW, FSPRx, M), col='blue')
#points(gW, M*( (1/(x^gW)-y)/(1-y) )^(1/gW), pch='.')
dev.off()


#SRR: over alphas and gammas to hit proxy
#SRR: over alphas and gammas to hit MSY
#SRR: add FSPRx & B40 visualizations
# one value of alpha to hit proxy (not MSY)
png(sprintf("srr%s_x%s_y%s.png", name, x, y))
#curve(SRR(x, getAlphaMsy(-1, FSPRx, M), getBeta(getAlphaMsy(-1, FSPRx, M), -1, M, 1), -1), lwd=3, xlim=c(0, 1), ylim=c(0, getAlphaMsy(-1, FSPRx, M)/getBeta(getAlphaMsy(-1, FSPRx, M), -1, M, 1)))
X=x
curve(SRR(x, getAlphaProxy(-1, M, X, y), getBeta(getAlphaProxy(-1, M, X, y), -1, M, 1), -1), 
	lty=2, 
	lwd=3,
	xlab="Biomass",
	ylab="Recruitment", #"Production",
	main=""
)
curve(SRR(x, getAlphaProxy(gW, M, X, y), getBeta(getAlphaProxy(gW, M, X, y), gW, M, 1), gW), add=T, lwd=3)
curve(M*x, col='red', add=T)
#curve(getAlphaMsy(-1, FSPRx, M)*x, add=T, lty=1)
#curve(getAlphaProxy(-1, M, X, y)*x, add=T, lty=2)
#curve(getAlphaProxy(gW, M, X, y)*x, add=T, lty=3)
segments(-1, M, 1, M, col='red')
segments(0, 0, 0.5, (M+FSPRx)*0.5, col='blue')
#
#FStar = FMsy(getAlphaMsy(-1, FSPRx, M), -1, M)
#BStar = PBar(getAlphaMsy(-1, FSPRx, M), getBeta(getAlphaMsy(-1, FSPRx, M), -1, M, 1), -1, FStar, M)
#abline(v=BStar, lty=1)
#
FStar = FMsy(getAlphaProxy(-1, M, X, y), -1, M)
BStar = PBar(getAlphaProxy(-1, M, X, y), getBeta(getAlphaProxy(-1, M, X, y), -1, M, 1), -1, FStar, M)
#segments(BStar, -1, BStar, (FStar+M)*BStar, lty=2)
rug(BStar, lwd=2, lty=2, ticksize=0.05)
#
FStar = FMsy(getAlphaProxy(gW, M, X, y), gW, M)
BStar = PBar(getAlphaProxy(gW, M, X, y), getBeta(getAlphaProxy(gW, M, X, y), gW, M, 1), gW, FStar, M)
#abline(v=BStar, lty=3)
segments(BStar, 0, BStar, (FStar+M)*BStar, col="blue")
rug(BStar, lwd=2, ticksize=0.05)
legend("bottomright", legend=c("BH at Proxy", "Schnute at Proxy at MSY"), lty=c(2, 1), lwd=3)
dev.off()

#
#M=0.1
#

#
M=0.1

#RP space: plot BH MSY with proxy point
#RP space: plot
i=1
y = ys[i]
x = xs[i]
#
FSPRx = M*((1/x)-1)
gW = getGammaProMsy(x, y)
#
png("rpProxAll0.1.png")
curve(1/(x+2), 0, 3, xlim=c(0, 3), ylim=c(0.2,0.55), lwd=3, xlab="F*/M", ylab="B*/B0", main="Reference Points")
points(FSPRx/M, y, col=i+1, pch=19)
#
FM = seq(0,3,length.out=100)
gLine = getZeta(gW, FM*M, M)
lines(FM, gLine, col=i+1, lwd=3)
#
gg = gW
for(i in 2:3){
	#
	x = xs[i]
	y = ys[i]
	#
	FSPRx = M*((1/x)-1)
	gW = getGammaProMsy(x, y)
	gg = c(gg, gW)
	#
	points(FSPRx/M, y, col=i+1, pch=19)
	gLine = getZeta(gW, FM*M, M)
	lines(FM, gLine, col=i+1, lwd=3)
}
legend("topright", legend=c(
	"BH MSY", 
	sprintf("%s Proxy (%s, %s)", names[1], xs[1], ys[1]), TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[1], round(gg[1],2))),
	sprintf("%s Proxy (%s, %s)", names[2], xs[2], ys[2]), TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[2], round(gg[2],2))), 
	sprintf("%s Proxy (%s, %s)", names[3], xs[3], ys[3]), TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[3], round(gg[3],2)))), 
	lwd=c(3, NA, 3, NA, 3, NA, 3), 
	pch=c(NA, 19, NA, 19, NA, 19, NA), 
	col=c(1, 2, 2, 3, 3, 4, 4)
)
dev.off()

#
#M=0.2
#

#
M=0.2

#RP space: plot BH MSY with proxy point
#RP space: plot
i=1
y = ys[i]
x = xs[i]
#
FSPRx = M*((1/x)-1)
gW = getGammaProMsy(x, y)
#
png("rpProxAll0.2.png")
curve(1/(x+2), 0, 3, xlim=c(0, 3), ylim=c(0.2,0.55), lwd=3, xlab="F*/M", ylab="B*/B0", main="Reference Points")
points(FSPRx/M, y, col=i+1, pch=19)
#
FM = seq(0,3,length.out=100)
gLine = getZeta(gW, FM*M, M)
lines(FM, gLine, col=i+1, lwd=3)
#
gg = gW
for(i in 2:3){
	#
	x = xs[i]
	y = ys[i]
	#
	FSPRx = M*((1/x)-1)
	gW = getGammaProMsy(x, y)
	gg = c(gg, gW)
	#
	points(FSPRx/M, y, col=i+1, pch=19)
	gLine = getZeta(gW, FM*M, M)
	lines(FM, gLine, col=i+1, lwd=3)
}
legend("topright", legend=c(
	"BH MSY", 
	sprintf("%s Proxy (%s, %s)", names[1], xs[1], ys[1]), TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[1], round(gg[1],2))),
	sprintf("%s Proxy (%s, %s)", names[2], xs[2], ys[2]), TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[2], round(gg[2],2))), 
	sprintf("%s Proxy (%s, %s)", names[3], xs[3], ys[3]), TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[3], round(gg[3],2)))), 
	lwd=c(3, NA, 3, NA, 3, NA, 3), 
	pch=c(NA, 19, NA, 19, NA, 19, NA), 
	col=c(1, 2, 2, 3, 3, 4, 4)
)
dev.off()

#
#table
#

M=0.2

#RP space: plot BH MSY with proxy point
#RP space: plot
i=1
y = ys[i]
x = xs[i]
Fy = c(0.5, 0.5, 2)
Bx = 1/((1/xs)-1+2) #c(0.29, 0.22, 0.21)
#
FSPRx = M*((1/x)-1)
gW = getGammaProMsy(x, y)
#
png("rpProxTable.png")
curve(1/(x+2), 0, 3, xlim=c(0, 3), ylim=c(0.2,0.55), lwd=3, xlab="F*/M", ylab="B*/B0", main="Reference Points")
points(FSPRx/M, y, col=i+1, pch=19)
segments(Fy[i], 1/(Fy[i]+2), FSPRx/M, y, col=i+1)
segments(FSPRx/M, Bx[i], FSPRx/M, y, col=i+1)
#
FM = seq(0,3,length.out=100)
gLine = getZeta(gW, FM*M, M)
#lines(FM, gLine, col=i+1, lwd=3)
#
gg = gW
for(i in 2:3){
	#
	x = xs[i]
	y = ys[i]
	#
	FSPRx = M*((1/x)-1)
	gW = getGammaProMsy(x, y)
	gg = c(gg, gW)
	#
	points(FSPRx/M, y, col=i+1, pch=19)
	segments(Fy[i], 1/(Fy[i]+2), FSPRx/M, y, col=i+1)
	segments(FSPRx/M, Bx[i], FSPRx/M, y, col=i+1)
	#gLine = getZeta(gW, FM*M, M)
	#lines(FM, gLine, col=i+1, lwd=3)
}
#
i=1
y=ys[i]
x=xs[i]
FSPRx = M*((1/x)-1)
segments(Fy[i], 1/(Fy[i]+2), FSPRx/M, y, col=i+1)
segments(FSPRx/M, Bx[i], FSPRx/M, y, col=i+1)
#
legend("topright", legend=c(
	"BH MSY", 
	sprintf("%s Proxy (%s, %s)", names[1], xs[1], ys[1]), #TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[1], round(gg[1],2))),
	sprintf("%s Proxy (%s, %s)", names[2], xs[2], ys[2]), #TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[2], round(gg[2],2))), 
	sprintf("%s Proxy (%s, %s)", names[3], xs[3], ys[3])), #TeX(sprintf("%s Schnute $\\gamma\\approx$%s", names[3], round(gg[3],2)))), 
	lwd=c(3, NA, NA, NA), 
	pch=c(NA, 19, 19, 19), 
	col=c(1, 2, 3, 4)
)
dev.off()









#alpha = 1
#gamma = 1#-1
#FMSY  = FMsy(alpha, gamma, M)
#FF = 0.4

##
#f = function(g){ M*( (1/(x^g)-y)/(1-y) )^(1/g) - getAlpha(g, FSPRx, M) }
#gN = uniroot(f,c(-1, 1))$root
##
#pow = y/(x-1)/(y-1)
#gW = pow - lambertWn(pow*x^pow*log(x))/log(x)

##proxy Alpha
#a = M*( (1/(x^gs)-y)/(1-y) )^(1/gs)
#aMsy = getAlpha(gs, FF, M)






#ugly = (y*x-y+1)/((1-x)*(1-y))
#gW = lambertWp( (log(x)*x^(ugly+1))/((1-x)*(1-y)) )/log(x) - ugly
#ugly = #(y*x-y+1)/((1-x)*(1-y))

##
#aP = uniroot(function(a){FProxyMsy(a, gamma, M, x)}, c(M,10))$root
#FaP = FMsy(aP, gamma, M)
#ZaP = getZeta(gamma, FaP, M)
##
#gP = uniroot(function(g){FProxyMsy(alpha, g, M, x)}, c(-2, 2))$root
#FgP = FMsy(alpha, gP, M)
#ZgP = getZeta(gP, alpha, M)

##
#FProxyMsy = function(alpha, gamma, M, x){
#	1/gamma - (1/gamma+1-x) * (M/(x*alpha))^gamma
#}

##
#alphas = seq(M+0.4, 10, length.out=100)
##
#gPs = c()
#FgPs = c()
#ZgPs = c()
#for(as in alphas){
#	#
#	gP = uniroot(function(g){FProxyMsy(as, g, M, x)}, c(-2, 2))$root
#	#
#	FgP = FMsy(as, gP, M)
#	ZgP = getZeta(gP, as, M)
#	#
#	FgPs = c(FgPs, FgP)
#	ZgPs = c(ZgPs, ZgP)
#	gPs = c(gPs, gP)
#}


#


##
##REGIMES
##
#
#cols = brewer.pal(9, 'Set1')
##cols = brewer.pal(9, 'Dark2')
#
##design
#png(sprintf("designLineColor%s.png", mod))
#plot(D$xiSeed, D$zetaSeed, xlab=TeX("$F_{MSY}/M$"), ylab=TeX('$B_{MSY}/B_0$'), main="Schnute RP Design", cex.lab = 1.5, cex.main= 1.5, xlim=c(xiBot, xiTop), ylim=c(zetaBot, zetaTop))
##curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#points(D$xiSeed[D$gSeed>=1], D$zetaSeed[D$gSeed>=1], col=cols[1], pch=19, cex=1.5)
#points(D$xiSeed[D$gSeed>=0 & D$gSeed<1], D$zetaSeed[D$gSeed>=0 & D$gSeed<1], col=cols[2], pch=19, cex=1.5)
#points(D$xiSeed[D$gSeed<0 & D$gSeed>=-1], D$zetaSeed[D$gSeed<0 & D$gSeed>=-1], col=cols[3], pch=19, cex=1.5)
#points(D$xiSeed[D$gSeed< -1], D$zetaSeed[D$gSeed< -1], col=cols[4], pch=19, cex=1.5)
#curve(1/(x+2), from=0, to=4, lwd=3, add=T)
##curve()
##legend('bottomleft', legend=c('Super-Logistic', 'Super-Ricker', 'Super-BH', 'Super-Cushing'), col=cols[1:4], pch=19)
#legend('bottomleft', legend=c('Logistic-Like', 'Ricker-Like', 'BH-Like', 'Cushing-Like'), col=cols[1:4], pch=19)
##dotX = lFXStar[!mask,2][freq]
##dotY = lFXStar[!mask,3][freq]
##points(dotX[dotY>zetaBot & dotY<zetaTop], dotY[dotY>zetaBot & dotY<zetaTop], pch='.')
###plot(D$xiSeed, D$zetaSeed, xlab=TeX("$F_{MSY}/M$"), ylab=TeX('$B_{MSY}/B_0$'), main="Design", cex.lab = 1.5, cex.main= 1.5)
###curve(1/(x+2), from=0, to=4, lwd=3, add=T)
#dev.off()

##
##DATA
##
#
##
#P0 = 10000
#mod = "HHardFlatT30N150WWideN112" #"HHardExpT45N150M0.05Wide" #"HHardFlatT30N150M0.1WWideN56" # #"HHardExpT45N150M0.3Wide" # #"HHardExpT45N150M0.1Wide" #"ExpT45N150M0.3Wide" #"HHardFlatT30N150WWideExtra" # #"ExpT45N150Wide" ##"ExpT30L4N150Wide" # #"ExpT45N150" #
#place = sprintf("./modsSchnute%s/", mod)
#
##
#xiRes = 0.5
#zetaTop = 0.6 #0.7
#zetaBot = 0.2 #0.1
#xiBot = 0.5
#xiTop = 3.5
#
##
#f = sprintf( "%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda"))[1] )
#M = readRDS(f)$M #0.2
#
##
#D = getData(place, c(xiBot, xiTop), c(zetaBot, 0.7))
#D = D[D$lFV>0 & D$lKV>0,]
#
##
#zetaStar = seq(min(D$zetaSeed), max(D$zetaSeed), 0.005) #length.out=)  #seq(0.15, 0.35, 0.001)  #rev(seq(0.1
#xiStar   = seq(min(D$xiSeed), max(D$xiSeed), 0.01) #length.out=)       #seq(1, 3.5, 0.005)
#lFXStar = cbind(1, expand.grid(xiStar, zetaStar))
#mask = sapply(1:nrow(lFXStar), function(i){
#               #
#               win = 0.3
#               #
#               xi = lFXStar[i,2]
#               zeta = lFXStar[i,3]
#               #step = win*10/0.1
#               bot = c()
#               top = c()
#               for(x in seq(max(xi-win, xiBot), min(xi+win, xiTop), length.out=30)){
#                       #
#                       isWin = x<=(D$xiBin+win) & x>=(D$xiBin-win)
#                       bot = c(bot, min(D$zetaBin[isWin]) )
#                       top = c(top, max(D$zetaBin[isWin]) )
#               }
#               bot = mean(bot)
#               top = mean(top)
#               ##
#               #bot = min(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#               #top = max(D$zetaSeed[D$xiSeed==round(xi/xiRes)*xiRes])
#               #
#               return( zeta>bot & zeta<top & zeta>zetaBot & zeta<zetaTop )
#               #return(T)
#       }
#)

