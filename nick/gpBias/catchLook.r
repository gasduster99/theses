rm(list=ls())

#
library(RColorBrewer)

#
TT = 31
tt = 1
time = tt:TT
mid = round(TT/2)
#
M = 0.2
cMax = 2 #4*M
cMin = 0.2 #M/4
slope = (cMax-cMin)/(mid-1)
#
b1 = cMin-slope*mid
b2 = 2*slope*mid + b1
cDown = (-slope*time+b2)*(time<=mid) + (slope*time+b1)*(time>mid)
#
b1 = cMax-slope*mid
b2 = 2*slope*mid + b1
cUp = (slope*time+b1)*(time<=mid) + (-slope*time+b2)*(time>mid)
#
b = log(cMin/cMax)/(tt-mid)
a = exp(log(cMax)-b*mid)
rSlope = (cMax-1)/(mid-1)
rb = 2*rSlope*mid + cMax-rSlope*mid
cExp = (a*exp(b*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
#
cFlat = rep(1, TT)
#
png("fishRates.png")
cols = brewer.pal(9, "Set1")
plot(time, cDown, 'l', lwd=3, col=cols[1], ylab="F(t)/F*")
lines(time, cUp, lwd=3, col=cols[2])
lines(time, cExp, lwd=2, col=cols[3])
lines(time, cFlat, lwd=3, col=cols[4])
dev.off()

