rm(list=ls())

#
library(latex2exp)

#
#FUNCTIONS
#

#
FMsy = function(alpha, gamma){ alpha/gamma }

#
#
#

##
#place = './contrastTestsX0.699Z0.201/'
#place = './contrastTestsX0.672Z0.323/'
#place = './contrastTestsX0.602Z0.214/'
#place = './contrastTestsX0.588Z0.157/'
#place = './contrastTestsX0.506Z0.246/'

#
j = 1
add = F
#
png('contrastTest3.png')
mins = c()
places = c('./contrastTestsX0.699Z0.201/', "./contrastTestsX0.244Z0.389/", "./contrastTestsX0.415Z0.506/") #"./contrastTestsX0.22Z0.413/") #, './contrastTestsX0.672Z0.323/', './contrastTestsX0.602Z0.214/', './contrastTestsX0.588Z0.157/', './contrastTestsX0.506Z0.246/')
for(place in places){
	#
	fitFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("fit*.rda")))
	datFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("dat*.rda")))
	#
	n = length(fitFiles)
	
	#
	D = data.frame(con=numeric(n), datFmsy=numeric(n), fitFmsy=numeric(n))
	for(i in 1:n){
		#
		dat = readRDS(datFiles[i])
		fit = readRDS(fitFiles[i])
		#
		D$con[i] = dat$con
		D$datFmsy[i] = FMsy(dat$alpha, dat$gamma)
		D$fitFmsy[i] = FMsy(fit$alpha, fit$gamma)
	}
	#rownames(D) = fitFiles
	
	#
	#plot(head(tail(D$con, 43), 15), head(tail(D$datFmsy-D$fitFmsy, 43), 15))
	who = D$fitFmsy<10 & D$con<=0.7 #0.8
	
	#
	if(add!=T){
		#
		xs = ( D$con[who]-min(D$con[who]) )/max(D$con[who])
		#png('contrastTest3.png')
		plot(xs, (D$datFmsy-D$fitFmsy)[who], 
			col=j, 
			ylim=c(-0.05,0.5),
			main=TeX("Bias in Estimated Schaefer $F_{MSY}$"),
			ylab="Bias",
			xlab="Contrast"
		)
		#dev.off()
		fs = c( dat$xi)
		bs = c( dat$zeta)
		mins = c( mins, D$con[who][(D$datFmsy-D$fitFmsy)[who]==min((D$datFmsy-D$fitFmsy)[who])] )
		add=T
	}else{
		fs = c(fs, dat$xi)
		bs = c(bs, dat$zeta)
		xs = ( D$con[who]-min(D$con[who]) )/max(D$con[who])
		points(xs, (D$datFmsy-D$fitFmsy)[who], col=c(1,2,4)[j])
		mins = c( mins, D$con[who][(D$datFmsy-D$fitFmsy)[who]==min((D$datFmsy-D$fitFmsy)[who])] )
	}
	#
	j=j+1
}
legend("topright", legend=c(TeX(sprintf("$F^*$=%.2f    $B^*/B_0$=%.2f", fs, bs))), col=c(1,2,4), pch=1)
dev.off()

