rm(list=ls())

#
#
#

#
hake  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
catch = c(94, 212, 195, 383, 320, 402, 366, 505, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
#					    606

#
png("hakeIndex.png")
plot(hake, xlab="Time", ylab="CPUE", main="Index Of Abundance")
dev.off()

#
png("hakeCatch.png")
plot(catch, xlab="Time", ylab="Biomass", main="Catch")
dev.off()

