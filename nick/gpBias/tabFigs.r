rm(list=ls())

#
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
#MAIN
#

##
#dir = "./modsPellaFineQFixReduxP010000" 
#M = 0.2
##
#dir = "./modsPellaFineQFixRFixP010000" 
#M = 0.2
#
mod = "CapNoQ"
dir = sprintf("./modsPella%s", mod)
M = 0.2

##
#catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)

#
xi = 1 		#2.3 	#
zeta = 0.2 	#0.25 	#
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
#grd = 0:fit$N0
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Recruitment",
	xlab="B",
	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(fit$N)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1, 
	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta)), 
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

#(1.0, 0.2)       &              &               &               &               &               &
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))


#
xi = 1
zeta = 0.6
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
#grd = 0:10000
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Recruitment",
	xlab="B",
	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(fit$N)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta)),
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
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))




#
xi = 3.5
zeta = 0.6
fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
#
dat = readRDS(fileDat)
fit = readRDS(fileFit)

#i
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
#grd = 0:fit$N0
lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))

#
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Recruitment",
	xlab="B",
	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(fit$N)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta)),
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
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))



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
yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)

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
png(sprintf("curveCompare%sX%sZ%s.png", mod, xi, zeta))
cols = brewer.pal(9, "Set1")
curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
	lwd=3, 
	ylim=c(0, yMax),
	ylab="Recruitment",
	xlab="B",
	main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta))
)
polygon(c(grd, rev(grd)), c(qLns[,1], rev(qLns[,3])), 
	col=adjustcolor(cols[1], alpha.f=0.2),
        border=F
)
curve(SRR(x, fit$lalpha, fit$lbeta, fit$gamma), 0, xMax, lwd=3, add=T, col=cols[1])
rug(fit$N)
legend("topright", legend=c("PT Truth", "Shaefer Fit"), col=c("black", cols[1]), lwd=3)
dev.off()

#
png(sprintf("fitCompare%sX%sZ%s.png", mod, xi, zeta))
plot(-1, -1,
        main=TeX(sprintf("$\\xi$=%s  $\\zeta$=%s", xi, zeta)),
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
ff = FMsy(fit$alpha, fit$gamma, M)
writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))


##
#xi = 3
#zeta = 0.5
#fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xi, zeta)
#fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xi, zeta)
##
#dat = readRDS(fileDat)
#fit = readRDS(fileFit)
#
##
#xMax = max(dat$N0, fit$N0)
#grd = 0:xMax
#yMax = max(SRR(grd, dat$lalpha, dat$lbeta, dat$gamma), SRR(grd, fit$lalpha, fit$lbeta, fit$gamma), na.rm=T)
#
##
#MM = 10^4
#who = c('lalpha', 'lbeta')
#C = fit$rsCov[who, who]
#m = c(fit$lalpha, fit$lbeta)
#sam = rmvnorm(MM, m, C)
##
#lns = mapply(function(a, b){SRR(grd, a, b, fit$gamma)}, sam[,1], sam[,2])
#qLns = rowQuantiles(lns, probs=c(0.025, 0.5, 0.975))
#
##
#png(sprintf("curveCompareX%sZ%s.png", xi, zeta))
#cols = brewer.pal(9, "Set1")
#curve(SRR(x, dat$lalpha, dat$lbeta, dat$gamma), 0, xMax, 
#	lwd=3, 
#	ylim=c(0, yMax),
#	ylab="Recruitment",
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
#
##
#ff = FMsy(fit$alpha, fit$gamma, M)
#writeLines(sprintf("(%1.1f, %1.1f)\t& %1.2f\t& %1.2f\t& %1.1f\t& %1.2f\t& %1.1e\t& %1.3f\t", dat$xi, dat$zeta, ff/M, PBar(ff, fit$alpha, fit$beta, fit$gamma, M)/fit$N0, fit$N0, fit$alpha, fit$q, fit$sdo))





