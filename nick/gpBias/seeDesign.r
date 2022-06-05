rm(list=ls())

#
#
#

#
place = './modsSchnuteFlatT30N150Wide/' #'./modsSchnuteExpN150/'

#
datFiles = sprintf("%s%s", place, list.files(path=place, pattern=glob2rx("datGen*.rda")))

#
xiList = c()
zetaList = c()
for(datF in datFiles){
	#
	dat = readRDS(datF)
	#
	print(c(dat$xi, dat$zeta))
	xiList = c(xiList, dat$xi)
	zetaList = c(zetaList, dat$zeta)
}

#
word = gsub("\\/", "", gsub("\\.", "", place))
png(sprintf('%sDesign.png', word))
plot(xiList, zetaList, pch=20)
#abline(v=xiL, lty=3)
#abline(h=zetaL, lty=3)
curve(1/(x+2), 0, 4, add=T)
dev.off()

