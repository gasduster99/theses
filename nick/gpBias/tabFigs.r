rm(list=ls())

#
library(pracma)
library(mvtnorm)
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
SRR = function(P, lalpha, lbeta, gamma){ 
	#
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
}

#
PBar = function(ff, alpha, beta, gamma, M){ beta*( 1 - (ff*(gamma-1)/alpha)^(1/(gamma-1)) ) }

#
FMsy = function(alpha, gamma, M){
        ##
        #FUpper = P0 #alpha-M
        ##root = uniroot.all(function(ff){  ((gamma-1)/alpha)^(gamma-1) - (M+ff)^(
        #root = uniroot.all(function(ff){ 1 - ((M+ff)*((gamma-1)/alpha))^(1/(gamma
        ##
        #if(length(root)<1){ return(NA) }
        #if(length(root)>1){ return(root[1]) }
        #
        root = alpha/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
        return(root)
}

#
diff = function(f, x, h=0.00001){
	(f(x)-f(x-h))/h
}

#
dPdt = function(t, P, lalpha, lbeta, gamma, M, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        Fmsy = exp(lalpha)/(gamma-1)*((gamma-1)/gamma)^(gamma-1)
        R = exp(lalpha)*P/(gamma-1)*(1-(P/exp(lbeta)))^(gamma-1)
        out = R - Fmsy*C*P
        #
        return( list(out) )
}

#
#MAIN
#

##
#dir = "./modsPellaFineQFixReduxP010000" 
#M = 0.2
##
#dir = "./modsPellaFineQFixRFixP010000" 
#M = 0.2
#
mod = "FlatNoQ"
dir = sprintf("./modsPella%s", mod)
M = 0.2
P0 = 10000

##
#catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)

#
xi = 0.5 		#2.3 	#
zeta = 0.2 	#0.25 	#
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
#grd = 0:fit$N0
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, xi, zeta))
diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
step=-10
quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
	step, step*diffDat,
	scale=100,
	col=1:length(dat$N) 
)
diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
	step, step*diffFit,
	scale=100,
	col=1:length(dat$N) 
)
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
dev.off()

###
##plot(diffDat, diffFit)
##lines(c(-100, 100), c(-100, 100))
##
#bioCol = brewer.pal(9, "YlGnBu")[c(-1, -2)]
##plot(diffDat[-1], diffFit[-1], col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
##lines(c(-100, 100), c(-100, 100))
#srrBar = mapply(function(a, b){SRR(fit$N, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSrr = rowQuantiles(srrBar, probs=c(0.025, 0.5, 0.975))
#plot(SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma), SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma), col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
#for(i in 1:nrow(qSrr)){
#	lines(rep(SRR(dat$N[i], dat$lalpha, dat$lbeta, dat$gamma), 2), c(qSrr[i,1], qSrr[i,3]), lwd=3, col=map2color(dat$N[i], bioCol, limits=c(0, 10000)))
#}
#lines(c(-10000, 10000), c(-10000, 10000))

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
NSam = iSam/fit$q
depSam = NSam/exp(sam[,2])
#depPSam = 

#
png(sprintf("bioCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Biomass",
        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N, lwd=3)
lines(fit$time, fit$N, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(NSam, probs=0.025),
			rev(colQuantiles(NSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
png(sprintf("depCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(depSam, probs=0.025),
			rev(colQuantiles(depSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
tmp = prodModel$new(
        dNdt=dPdt, N0Funk=function(lbeta){exp(lbeta)}, #model
        time=fit$time, catch=fit$catch, M=M,                          #con
        alpha=fit$alpha, beta=P0, gamma=2, #parameters
        lalpha=log(fit$alpha), lbeta=log(P0),       #reparameterize
        lq=log(0.00049), lsdo=log(0.01160256) #log(0.1160256)
)
dPost = t(mapply( function(a, b){
	#
	tmp$lalpha = a
	tmp$lbeta  = b
	tmp$iterate('lsode')
	#
	return( tmp$N/tmp$N[1] )
}, sam[,1], sam[,2] ))

#
png(sprintf("depPostCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, colMeans(dPost), lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(dPost, probs=0.025),
			rev(colQuantiles(dPost, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()


#(1.0, 0.2)       &              &               &               &               &               &
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
bl = ff*fit$catch*fit$N
dbl = dff*dat$catch*dat$N
lbl = TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))

#
xi = 0.5
zeta = 0.6
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
#grd = 0:10000
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, xi, zeta))
diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
step=-10
quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
	step, step*diffDat,
	scale=100,
	col=1:length(dat$N) 
)
diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
	step, step*diffFit,
	scale=100,
	col=1:length(dat$N) 
)
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
dev.off()

##
#plot(diffDat[-1], diffFit[-1])
#lines(c(-100, 100), c(-100, 100))

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
NSam = iSam/fit$q
depSam = NSam/exp(sam[,2])

#
png(sprintf("bioCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Biomass",
        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N, lwd=3)
lines(fit$time, fit$N, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(NSam, probs=0.025),
			rev(colQuantiles(NSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
png(sprintf("depCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(depSam, probs=0.025),
			rev(colQuantiles(depSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

dPost = t(mapply( function(a, b){
	#
	tmp$lalpha = a
	tmp$lbeta  = b
	tmp$iterate('lsode')
	#
	return( tmp$N/tmp$N[1] )
}, sam[,1], sam[,2] ))

#
png(sprintf("depPostCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, colMeans(dPost), lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(dPost, probs=0.025),
			rev(colQuantiles(dPost, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()


#
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
tl = ff*fit$catch*fit$N
dtl = dff*dat$catch*dat$N
ltl = TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))


#
xi = 3.5
zeta = 0.6
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
#grd = 0:fit$N0
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, xi, zeta))
diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
step=-10
quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
	step, step*diffDat,
	scale=100,
	col=1:length(dat$N) 
)
diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
	step, step*diffFit,
	scale=100,
	col=1:length(dat$N) 
)
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
dev.off()

##
#bioCol = brewer.pal(9, "YlGnBu")[c(-1, -2)]
##plot(diffDat[-1], diffFit[-1], col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
##lines(c(-100, 100), c(-100, 100))
#plot(SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma), SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma), col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
#lines(c(-10000, 10000), c(-10000, 10000))
##
#bioCol = brewer.pal(9, "YlGnBu")[c(-1, -2)]
##plot(diffDat[-1], diffFit[-1], col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
##lines(c(-100, 100), c(-100, 100))
#srrBar = mapply(function(a, b){SRR(fit$N, a, b, fit$gamma)}, sam[,1], sam[,2])
#qSrr = rowQuantiles(srrBar, probs=c(0.025, 0.5, 0.975))
#plot(SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma), SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma), col=map2color(dat$N, bioCol, limits=c(0, 10000)), pch=19)
#for(i in 1:nrow(qSrr)){
#	lines(rep(SRR(dat$N[i], dat$lalpha, dat$lbeta, dat$gamma), 2), c(qSrr[i,1], qSrr[i,3]), lwd=3, col=map2color(dat$N[i], bioCol, limits=c(0, 10000)))
#}
#lines(c(-10000, 10000), c(-10000, 10000))




#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
NSam = iSam/fit$q
depSam = NSam/exp(sam[,2])

#
png(sprintf("bioCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Biomass",
        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N, lwd=3)
lines(fit$time, fit$N, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(NSam, probs=0.025),
			rev(colQuantiles(NSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
png(sprintf("depCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(depSam, probs=0.025),
			rev(colQuantiles(depSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

dPost = t(mapply( function(a, b){
	#
	tmp$lalpha = a
	tmp$lbeta  = b
	tmp$iterate('lsode')
	#
	return( tmp$N/tmp$N[1] )
}, sam[,1], sam[,2] ))

#
png(sprintf("depPostCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, colMeans(dPost), lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(dPost, probs=0.025),
			rev(colQuantiles(dPost, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()


#
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
tr = ff*fit$catch*fit$N
dtr = dff*dat$catch*dat$N
ltr = TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))

#
xi = 3.5
zeta = 0.2
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
#grd = 0:xMax
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), qLns, na.rm=T)

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, xi, zeta))
diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
step=-10
quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
	step, step*diffDat,
	scale=100,
	col=1:length(dat$N) 
)
diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
	step, step*diffFit,
	scale=100,
	col=1:length(dat$N) 
)
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
iSam = sapply(fit$N, function(x){rlnorm(MM*100, fit$lq+log(x), exp(fit$lsdo))})
NSam = iSam/fit$q
depSam = NSam/exp(sam[,2])

#
png(sprintf("bioCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Biomass",
        ylim=c(min(0, dat$N), max(fit$N, dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N, lwd=3)
lines(fit$time, fit$N, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(NSam, probs=0.025),
			rev(colQuantiles(NSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
png(sprintf("depCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, fit$N/fit$N0, lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(depSam, probs=0.025),
			rev(colQuantiles(depSam, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

dPost = t(mapply( function(a, b){
	#
	tmp$lalpha = a
	tmp$lbeta  = b
	tmp$iterate('lsode')
	#
	return( tmp$N/tmp$N[1] )
}, sam[,1], sam[,2] ))

#
png(sprintf("depPostCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="Depletion",
        ylim=c(min(0, dat$N/dat$N0), max(fit$N/fit$N0, dat$N/dat$N0, colQuantiles(depSam, probs=0.975))),
        xlim=c(min(dat$time), max(dat$time))
)
lines(dat$time, dat$N/dat$N0, lwd=3)
lines(fit$time, colMeans(dPost), lwd=3, col='red')
polygon( c(fit$time, rev(fit$time)),
                c(
			colQuantiles(dPost, probs=0.025),
			rev(colQuantiles(dPost, probs=0.975))
		),
                col = makeTransparent('red', alpha=100),
                border = NA
)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()


#
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))
br = ff*fit$catch*fit$N
dbr = dff*dat$catch*dat$N
lbr = TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))

#
xi = 3.5
zeta = 0.5
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

##
#png(sprintf("curveCompareArrow%sX%sZ%s.png", mod, xi, zeta))
#diffDat = diff(function(x){SRR(x, dat$lalpha, dat$lbeta, dat$gamma)}, dat$N)
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
#)
#step=-10
#quiver( dat$N, SRR(dat$N, dat$lalpha, dat$lbeta, dat$gamma),
#	step, step*diffDat,
#	scale=100,
#	col=1:length(dat$N) 
#)
#diffFit = diff(function(x){SRR(x, fit$lalpha, fit$lbeta, fit$gamma)}, fit$N)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#quiver( fit$N, SRR(fit$N, fit$lalpha, fit$lbeta, fit$gamma),
#	step, step*diffFit,
#	scale=100,
#	col=1:length(dat$N) 
#)
#dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))



#
xi = 2
zeta = 0.5
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))



#
xi = 1
zeta = 0.5
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#
xMax = max(dat$N0, fit$N0)
grd = 0:xMax
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

#
MM = 10^4
who = c('lalpha', 'lbeta')
C = fit$rsCov[who, who]
m = c(fit$lalpha, fit$lbeta)
sam = rmvnorm(MM, m, C)
#
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, n=1000,
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Production",
	xlab="B",
	main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(dat$N, line=0.75)
rug(fit$N, col=cols[1])
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

#
dff = FMsy(dat$alpha, dat$gamma, M)
bMsy = PBar(dff, dat$alpha, dat$beta, dat$gamma, M)
#
png(sprintf("fitCompareBmsy%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$F^*$=%s  \t  $B^*/B_0$=%s", xi*M, zeta)),
        xlab="Time",
        ylab="CPUE",
        ylim=c(min(0, dat$q*dat$N), max(fit$q*fit$N, dat$q*dat$N)),
        xlim=c(min(dat$time), max(dat$time))
)
dat$plotMean(add=T)
fit$plotMean(add=T, col='red')
fit$plotBand(col='red')
abline(h=bMsy*dat$q, lty=2)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", 'red'), lwd=3)
dev.off()

##
#png(sprintf("catch%s.png", mod))
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, 1, 0.5)
#dat = readRDS(fileDat)
#plot(dat$time, dat$catch, xlab="Time", ylab="F(t)/F*", ylim=c(0, 2), main="Catch=F* (F(t)/F*) B(t)", type="l", lwd=3)
#dev.off()

#
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))

#
png(sprintf("fts%s.png", mod))
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, 1, 0.5)
dat = readRDS(fileDat)
#plot(dat$time, dat$catch, xlab="Time", ylab="F(t)/F*", ylim=c(0, 2), main="Catch=F* (F(t)/F*) B(t)", type="l", lwd=3)
plot(dat$time, dbl, type='l', col=cols[1], lwd=3, xlab="Time", ylab="Catch", ylim=c(min(bl, tl, tr, br), max(bl, tl, tr, br)), main="Catch")
lines(dat$time, bl, col=cols[1], lty=2, lwd=3)
lines(dat$time, dtl, col=cols[2], lwd=3)
lines(dat$time, tl, col=cols[2], lty=2, lwd=3)
lines(dat$time, dtr, col=cols[3], lwd=3)
lines(dat$time, tr, col=cols[3], lty=2, lwd=3)
lines(dat$time, dbr, col=cols[4], lwd=3)
lines(dat$time, br, col=cols[4], lty=2, lwd=3)
legend("topright", legend=c(lbl, ltl, ltr, lbr), col=cols, lwd=3)
dev.off()








##
#png(sprintf("curveCompareX%sZ%s.png", xi, zeta))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Production",
#	xlab="B",
#	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta))
#)
#polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
#	col=adjustcolor(cols[1], alpha.f=0.2),
#        border=F
#)
#curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
#legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
#dev.off()

##
#ff = FMsy(fit$alpha, fit$gamma, M)
#writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))





