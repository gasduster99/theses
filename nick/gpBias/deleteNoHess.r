rm(list=ls())

#
library(VGAM)
library(pracma)
library(mvtnorm)
library(plot.matrix)
library(RColorBrewer)
library(scatterplot3d)

#
#FUNCTIONS
#

#
dist = function(x, xi, zeta){ sqrt((xi-x)^2 + (zeta-1/(x+2))^2) }

#
getData = function(dir, xiSims, zetaSims){
        #dir    : a directory containing data
        #xiSims : xis of simulated data
        #zetaSims: zetas of simulated data

        #
        di = 1
        D = data.frame(xi=double(), zeta=double(), xiHat=double(), xiMin=double(), zetaHat=double(), distHat=double(), distMin=double(), cllhs=double(), cllhsV=double(), stringsAsFactors=F)
        for(i in 1:length(zetaSims)){
                for(j in 1:length(xiSims)){
                        #
                        fileDat = sprintf('%s/datGen_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])
                        fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])

                        #
                        if( file.exists(fileDat) & file.exists(fileFit) ){
                                #
                                dat = readRDS(fileDat)
                                fit = readRDS(fileFit)

                                #
                                xiHat   = -log(1-cloglog(fit$cllhs, inverse=T))/M
                                zetaHat = 1/(xiHat+2)
                                distHat = sqrt((dat$xi-xiHat)^2 + (dat$zeta-zetaHat)^2)
                                opt = optimize(dist, c(0, dat$xi), xi=dat$xi, zeta=dat$zeta)
                                #
                                if(length(fit$rsCov)==0){ v=0 }else{ v=fit$rsCov['cllhs', 'cllhs'] }
                                D[di,]  = c(dat$xi, dat$zeta, xiHat, opt$minimum, zetaHat, distHat, opt$objective, fit$cllhs, v)
                                di=di+1
                        }
                }
        }
        #
        return(D)
}

rmNoHess = function(dir, xiSims, zetaSims){
        #dir    : a directory containing data
        #xiSims : xis of simulated data
        #zetaSims: zetas of simulated data

	#	
        for(i in 1:length(zetaSims)){
                for(j in 1:length(xiSims)){
                        #
                        fileFit = sprintf('%s/fit_xi%s_zeta%s.rda', dir, xiSims[j], zetaSims[i])

                        #
                        if( file.exists(fileFit) ){
                                #
                                fit = readRDS(fileFit)
				
                                #
                                if(length(fit$rsCov)==0){ system(sprintf('rm %s', fileFit)) }
			}
		}
	}
}

#
#DATA
#

#
dir = './modsFine/'
#
zetaSims = rev(seq(0.1, 0.75, 0.01)) #rev(seq(0.25, 0.75, 0.01))
xiSims =   seq(0.5, 3.5, 0.05)
#
rmNoHess(dir, xiSims, zetaSims)



