rm(list=ls())

#
library(geoR)
library(boot)
library(scatterplot3d)

#
#FUNCTIONS
#

#
boxTrans = function(x, lambda, lambda2=0){
        if(lambda==0){return(log(x+lambda2))}
        return( ((x+lambda2)^lambda-1)/lambda )
}

#
dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/(x+2))^2) }

#
K = function(x1, x2){ exp((-(x1-x2)^2)/(2*l^2)) }

#
#MAIN
#

#
M = 0.2
#
zetaSims = rev(seq(0.25, 0.75, 0.1))
xiSims = seq(0.6, 3.1, 0.5)

#
i = 1
D = data.frame(xi=double(), xiHat=double(), xiMin=double(), zeta=double(), zetaHat=double(), distHat=double(), distMin=double(), stringsAsFactors=F)
for(zi in 1:length(zetaSims)){
        for(xj in 1:length(xiSims)){
		#
		fileDat = sprintf('mods/datGen_xi%s_zeta%s.rda', xiSims[xj], zetaSims[zi])
		fileMod = sprintf('mods/bevH_xi%s_zeta%s.rda', xiSims[xj], zetaSims[zi])

		if( file.exists(fileDat) & file.exists(fileMod) ){
			#
			dat = readRDS(fileDat)
			mod = readRDS(fileMod)

			#
			xiHat   = -log(1-inv.logit(mod$lhs))/M
			zetaHat = 1/(mod$xi+2)
			distHat = sqrt((dat$xi-xiHat)^2 + (dat$zeta-zetaHat)^2)
			opt = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)
			#
			D[i,]  = c(dat$xi, xiHat, opt$minimum, dat$zeta, zetaHat, distHat, opt$objective)
			i=i+1
		}	
	}
}

#
pdf('distanceTransSmall.pdf', height=15, width=10)
layout(matrix(1:6, nrow=3, ncol=2, byrow=T))
cex = 1
#
p = scatterplot3d(D$xi, D$zeta, D$distHat, cex.lab=cex, cex.axis=cex)
plane = lm(distHat~xi+zeta, data=D)
p$plane3d(plane)
#
pp = predict(plane)
y = D$distHat-pp
scatterplot3d(D$xi, D$zeta, D$distHat-pp)

#
p = scatterplot3d(D$xi, D$zeta, D$distHat-D$distMin)
plane = lm(distHat-distMin~xi+zeta, data=D)
p$plane3d(plane)
#
pp = predict(plane)
y = (D$distHat-D$distMin)-pp
scatterplot3d(D$xi, D$zeta, (D$distHat-D$distMin)-pp)


#
p = scatterplot3d(D$xi, D$zeta, log(D$distHat-D$distMin))
plane = lm(log(distHat-distMin)~xi+zeta, data=D)
p$plane3d(plane)
#
pp = predict(plane)
y = log(D$distHat-D$distMin)-pp
scatterplot3d(D$xi, D$zeta, log(D$distHat-D$distMin)-pp)

#
dev.off()


