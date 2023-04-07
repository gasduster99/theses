rm(list=ls())

#
library(pracma)
library(mvtnorm)
library(rootSolve)
library(latex2exp)
library(matrixStats)
library(RColorBrewer)
#
source('prodClass0.1.1.r')

#
#FUNCTIONS
#

#
map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#
#pellaSRR = function(P, alpha, beta, gamma){ alpha*P/(gamma-1)*(1-(P/beta)^(gamma-1)) }
##
#SRR = function(B, alpha, beta, gamma){
#        alpha*B*(1-beta*gamma*B)^(1/gamma)
#}
SRR = function(B, lalpha, lbeta, gamma){
        ##Schnute
	#exp(lalpha)*B*(1-exp(lbeta)*gamma*B)^(1/gamma)
	#PT
	exp(lalpha)*B/(gamma-1)*(1-(B/exp(lbeta))^(gamma-1))
}

#
#PLOT
#

#
EM = 10000
#
mod = "ExpT45" #"FlatT30" #
dir = sprintf("./modsPT%s/", mod)
P0 = 10000

#LHS boundaries
xiMin = 0.1
xiMax = 0.7 #2 #3.5
zetaMin = 0.2
zetaMax = 0.6
#
datFiles = sprintf("%s%s", dir, list.files(path=dir, pattern=glob2rx("datGen*.rda")))
#
not = c(
	"./modsPTExpT45/datGen_xi0.079_zeta0.187.rda",
	"./modsPTExpT45/datGen_xi0.094_zeta0.202.rda",
	"./modsPTExpT45/datGen_xi0.137_zeta0.218.rda",
	"./modsPTExpT45/datGen_xi0.119_zeta0.225.rda",
	"./modsPTExpT45/datGen_xi0.702_zeta0.614.rda",
	"./modsPTExpT45/datGen_xi0.741_zeta0.592.rda",
	"./modsPTExpT45/datGen_xi0.735_zeta0.573.rda",
	"./modsPTExpT45/datGen_xi0.75_zeta0.597.rda",
	"./modsPTExpT45/datGen_xi0.649_zeta0.612.rda",
	"./modsPTExpT45/datGen_xi0.682_zeta0.546.rda",
	"./modsPTExpT45/datGen_xi0.696_zeta0.533.rda",
	"./modsPTExpT45/datGen_xi0.634_zeta0.583.rda",
	"./modsPTExpT45/datGen_xi0.77_zeta0.61.rda",
	"./modsPTExpT45/datGen_xi0.774_zeta0.586.rda",
	"./modsPTExpT45/datGen_xi0.706_zeta0.678.rda",
	"./modsPTExpT45/datGen_xi0.758_zeta0.656.rda",
	"./modsPTExpT45/datGen_xi0.718_zeta0.679.rda",
	"./modsPTExpT45/datGen_xi0.781_zeta0.594.rda",
	"./modsPTExpT45/datGen_xi0.675_zeta0.686.rda",
	"./modsPTExpT45/datGen_xi0.712_zeta0.689.rda",
	"./modsPTExpT45/datGen_xi0.788_zeta0.636.rda",
	"./modsPTExpT45/datGen_xi0.778_zeta0.545.rda",
	"./modsPTExpT45/datGen_xi0.612_zeta0.535.rda",
	"./modsPTExpT45/datGen_xi0.787_zeta0.672.rda",
	"./modsPTExpT45/datGen_xi0.611_zeta0.517.rda", #
	"./modsPTExpT45/datGen_xi0.663_zeta0.483.rda",
	"./modsPTExpT45/datGen_xi0.657_zeta0.477.rda",
	"./modsPTExpT45/datGen_xi0.58_zeta0.658.rda",
	"./modsPTExpT45/datGen_xi0.705_zeta0.463.rda",
	"./modsPTExpT45/datGen_xi0.597_zeta0.507.rda",
	"./modsPTExpT45/datGen_xi0.568_zeta0.66.rda",
	"./modsPTExpT45/datGen_xi0.56_zeta0.64.rda",
	"./modsPTFlatT30/datGen_xi0.678_zeta0.612.rda",	
	"./modsPTFlatT30/datGen_xi0.103_zeta0.586.rda",
	"./modsPTFlatT30/datGen_xi0.085_zeta0.611.rda",
	"./modsPTFlatT30/datGen_xi0.079_zeta0.616.rda",
	"./modsPTFlatT30/datGen_xi0.082_zeta0.623.rda"
)
datFiles = datFiles[!(datFiles%in%not)]

#
n = length(datFiles)
minDiff = min((zetaMax-zetaMin)/n, (xiMax-xiMin)/n)
binTrk  = ceiling(abs(log10(minDiff)))

#
xis = unlist(sapply(datFiles, function(fn){
                dat = readRDS(fn)
                if( length(dat$zeta)!=0 & length(dat$xi)!=0){
                        return(dat$xi)
                }
        })
)
zetas = unlist(sapply(datFiles, function(fn){
                dat = readRDS(fn)
                if( length(dat$zeta)!=0 & length(dat$xi)!=0){
                        return(dat$zeta)
                }
        })
)

###############################################################################################
##
##MAX, MAX
##
###############################################################################################

#minimize (max, max) norm
norms = sqrt((xiMax-xis)^2 + (zetaMax-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]

#
fileDat = names(who)
fileFit = gsub("datGen", "fit", fileDat)

#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = seq(0, xMax, length.out=EM)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, yMax),
        ylab="Production",
        xlab="B",
        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk), round(zeta,binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, lwd=3, ticksize=0.05)
#rug(dat$N, line=0.75)
#rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

###############################################################################################
##
##MIN, MAX
##
###############################################################################################

#minimize (max, max) norm
norms = sqrt((xiMin-xis)^2 + (zetaMax-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]

#
fileDat = names(who)
fileFit = gsub("datGen", "fit", fileDat)

#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = seq(0, xMax, length.out=EM)
#grd = 0:xMax

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, yMax),
        ylab="Production",
        xlab="B",
        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk), round(zeta,binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, lwd=3, ticksize=0.05)
#rug(dat$N, line=0.75)
#rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

###############################################################################################
##
##MIN, MIN
##
###############################################################################################

#minimize (max, max) norm
norms = sqrt((xiMin-xis)^2 + (zetaMin-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]

#
fileDat = names(who)
fileFit = gsub("datGen", "fit", fileDat)

#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = seq(0, xMax, length.out=EM)
#grd = 0:xMax

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, yMax),
        ylab="Production",
        xlab="B",
        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk), round(zeta,binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, lwd=3, ticksize=0.05)
#rug(dat$N, line=0.75)
#rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

###############################################################################################
##
##MAX, Min
##
###############################################################################################

#minimize (max, max) norm
norms = sqrt((xiMax-xis)^2 + (zetaMin-zetas)^2)
who   = which(min(norms)==norms)
xi    = xis[who]
zeta  = zetas[who]

#
fileDat = names(who)
fileFit = gsub("datGen", "fit", fileDat)

#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = seq(0, xMax, length.out=EM)
#grd = 0:xMax

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qSRR = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#qYield = rowQuantiles(lns-M*grd, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
png(sprintf("srrCompare%sX%sZ%s.png", mod, round(xi, binTrk), round(zeta, binTrk)))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
        lwd=3,
        ylim=c(0, yMax),
        ylab="Production",
        xlab="B",
        main=TeX(sprintf("$F_{MSY}$=%s  \t  $B_{MSY}/B_0$=%s", round(xi, binTrk), round(zeta,binTrk))),
        cex.lab = 1.5,
        cex.main= 1.5
)
polygon(c(grd, rev(grd)), c(qSRR[,1], rev(qSRR[,3])),
        col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, lwd=3, ticksize=0.05)
#rug(dat$N, line=0.75)
#rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "BH Fit"), col=c("black", cols[1]), lwd=3)
dev.off()


