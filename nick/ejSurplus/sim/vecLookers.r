rm(list=ls())

#
library(TTR)
library(RColorBrewer)

#
#FUNCTIONS
#

#
addalpha <- function(colors, alpha=1.0) {
	r <- col2rgb(colors, alpha=T)
	# Apply alpha
	r[4,] <- alpha*255
	r <- r/255.0
	return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

#
#RESTRUCTURE
#

#
files = Sys.glob('bRat_*_M_0.04_0.30DatSD.csv')
#
MTrue = seq(0.05, 0.215, 0.007)
aList = list()

#Ms = list()
#yFills = list()
#xFills = list()
##
for(file in files){
	#
	dat = read.csv(file)
	dat = dat[order(dat$M),]
	#
	head.x = dat$median.x #dat$mean.x 
	head.y = dat$median.y #dat$mean.y 	
	tail.x = dat$tail.x
	tail.y = dat$tail.y
	Ms = dat$M
	#
	aList[[file]] = matrix(NA, nrow=length(MTrue), ncol=5)
	colnames(aList[[file]]) = c('M', 'tail.x', 'tail.y', 'head.x', 'head.y')
	#
	aList[[file]][,1] = MTrue
	for(M in Ms){
		#
		where = Ms==M
		#
		aList[[file]][where,'tail.x'] = tail.x[where]
		aList[[file]][where,'tail.y'] = tail.y[where]
		aList[[file]][where,'head.x'] = head.x[where]
		aList[[file]][where,'head.y'] = head.y[where]
	}
}

#
#PLOT
#

#
dRatBH = seq(0.05, 20, 0.001)
#for(g in glob){ dRatBH = c(dRatBH, as.numeric(strsplit(g, "_")[[1]][2])) }
bRatBH = dRatBH/(dRatBH*(dRatBH+2))
#
bRats = seq(0.3, 0.6, 0.1)
dRats = seq(0.5, 3, 0.5)
xcols = brewer.pal(length(bRats), 'Spectral')
ycols = brewer.pal(length(dRats), 'Spectral')
for(m in 1:length(MTrue)){
	#
	png(sprintf("./pics/full_M_0.04_0.30Gif%1.2fSD.png", MTrue[m]))
	plot(dRatBH, bRatBH, 'l', lwd=3, 
		xlim=c(0, 20), 
	        xaxs="i",
	        ylim=c(0,1),
		main=sprintf("M=%1.2f", MTrue[m]),
                xlab=expression(F[msy]/M),
                ylab=expression(B[msy]/B[0])
	)
	#i = 1
	for(file in files){	
		#
		pars = as.numeric(strsplit(file, "_")[[1]][2])
		xwhere = which(pars==bRats)
		arrows(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'], x1=aList[[file]][m,'head.x'], y1=aList[[file]][m,'head.y'],
		        col=addalpha(xcols[xwhere], 0.7),
		        lwd=1,
		        length=0.15,
		        xlim=c(0, 20) 
		        #xaxs="i",
		        #ylim=c(0,0.5)
		)
		points(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'],
		        pch=19,
		        col=addalpha(xcols[xwhere], 0.7)
		)
		#i = i+1
	}
	for(file in files){	
		#
		pars = as.numeric(strsplit(file, "_")[[1]][4])
                ywhere = which(pars==dRats)
		arrows(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'], x1=aList[[file]][m,'head.x'], y1=aList[[file]][m,'head.y'],
		        col=addalpha(ycols[ywhere], 0.5),
		        lwd=1,
		        length=0.15,
		        xlim=c(0, 20) 
		        #xaxs="i",
		        #ylim=c(0,0.5)
		)
		points(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'],
		        pch=19,
		        col=addalpha(ycols[ywhere], 0.5)
		)
		#i = i+1
	}
	dev.off()
}
#
system(sprintf('convert -delay 40 ./pics/full_M_0.04_0.30Gif*SD.png ./pics/full_M_0.04_0.30SD.gif'))

#
pdf(sprintf("./full_M_0.04_0.30GridSD.pdf"))
layout(matrix(1:12, 3, 4, byrow=T))
i = 1
for(m in rev(c(1, 8, length(MTrue)))){
	#
	j = 1
	for(b in bRats){
		#
		par(mar=c(2, 2, 3, 0.5))
                plot(-1, -1, 'l', lwd=3,
                        xlim=c(0, 20),
                        xaxs="i",
                        ylim=c(0,1),
                        #main=expression(B[msy]/B[0]*sprintf(" = %1.1f", b)), #sprintf("M=%1.2f", MTrue[m]),
                        xlab=expression(F[msy]/M),
                        ylab=expression(B[msy]/B[0])
                )
		#
		for(file in Sys.glob(sprintf('bRat_%1.1f*_M_0.04_0.30DatSD.csv', b))){
        	        ##
        	        #pars = as.numeric(strsplit(file, "_")[[1]][2])
        	        #xwhere = which(pars==bRats)
        	        arrows(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'], x1=aList[[file]][m,'head.x'], y1=aList[[file]][m,'head.y'],
        	                col=addalpha(xcols[j], 0.7),#xcols[i], 
        	                lwd=1,
        	                length=0.15,
        	                xlim=c(0, 20)
        	                #xaxs="i",
        	                #ylim=c(0,0.5)
        	        )
        	        points(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'],
        	                pch=20,
        	                col=addalpha(xcols[j], 0.7)  #xcols[i], 
        	        )
        	        #i = i+1
        	}
		for(file in Sys.glob(sprintf('bRat_%1.1f*_M_0.04_0.30DatSD.csv', b))){
        	        ##
        	        #pars = as.numeric(strsplit(file, "_")[[1]][2])
        	        #xwhere = which(pars==bRats)
        	        arrows(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'], x1=aList[[file]][m,'head.x'], y1=aList[[file]][m,'head.y'],
        	                col=addalpha(xcols[i], 0.7), #xcols[i], 
        	                lwd=1,
        	                length=0.15,
        	                xlim=c(0, 20)
        	                #xaxs="i",
        	                #ylim=c(0,0.5)
        	        )
        	        points(aList[[file]][m,'tail.x'], aList[[file]][m,'tail.y'],
        	                pch=20,
        	                col=addalpha(xcols[i], 0.7)  #xcols[i]
        	        )
        	        #i = i+1
        	}
		j = j+1
	}
	i = i+1
}
dev.off()





#
#BONE YARD
#

##effect first Column
		#if(b==bRats[1]){
		#	#
		#	par(mar=c(5, 3, 1, 0.5))
        	#	plot(-1, -1, 'l', lwd=3,
        	#	        xlim=c(0, 20),
        	#	        xaxs="i",
        	#	        ylim=c(0,1),
        	#	        main=expression(B[msy]/B[0]*sprintf(" = %1.1f", b)), #sprintf("M=%1.2f", MTrue[m]),
        	#	        xlab=expression(F[msy]/M),
        	#	        ylab=expression(B[msy]/B[0])
        	#	)
		##mtext(expression(B[msy]/B[0]*sprintf(" = %1.1f", b)), side=1, line=4)
		#}else if(b==bRats[3]){
		#	#
		#	par(mar=c(5, 3, 0.5, 0.5))
        	#	plot(-1, -1, 'l', lwd=3,
        	#	        xlim=c(0, 20),
        	#	        xaxs="i",
        	#	        ylim=c(0,1),
        	#	        main=expression(B[msy]/B[0]*sprintf(" = %1.1f", b)), #sprintf("M=%1.2f", MTrue[m]),
        	#	        xlab=expression(F[msy]/M),
        	#	        ylab=expression(B[msy]/B[0])
        	#	)
		#}else{
		#	#
		#	par(mar=c(5, 3, 0.5, 0.5))
        	#	plot(-1, -1, 'l', lwd=3,
        	#	        xlim=c(0, 20),
        	#	        xaxs="i",
        	#	        ylim=c(0,1),
        	#	        main=expression(B[msy]/B[0]*sprintf(" = %1.1f", b)), #sprintf("M=%1.2f", MTrue[m]),
        	#	        xlab=expression(F[msy]/M),
        	#	        ylab=expression(B[msy]/B[0])
        	#	)
		#}
#	#
#	smaX5 = SMA(head.x-dat$tail.x, n=5)
#	fillX = c(
#		(head.x-dat$tail.x)[1],
#		SMA(head.x-dat$tail.x, n=2)[2],
#	        SMA(head.x-dat$tail.x, n=3)[3],
#	        SMA(head.x-dat$tail.x, n=4)[4],
#		na.omit(smaX5) #na.omit(smaX5)[1]
#	)
#	#
#	smaY5 = SMA(head.y-dat$tail.y, n=5)
#	fillY = c(
#	        (head.y-dat$tail.y)[1],
#	        SMA(head.y-dat$tail.y, n=2)[2],
#	        SMA(head.y-dat$tail.y, n=3)[3],
#	        SMA(head.y-dat$tail.y, n=4)[4],
#	        na.omit(smaY5)#na.omit(smaY5)[1]
#	)
#	#
#	Ms[[file]] = dat$M
#	xFills[[file]] = fillX
#	yFills[[file]] = fillY
#}
##
#cols = rep(brewer.pal(11, 'Spectral')[-6], 2)
##
#dev.new()
#plot(Ms[[1]], yFills[[1]], 'l', 
#	ylim=c(do.call(min, yFills), do.call(max, yFills)),
#	xlim=c(do.call(min, Ms), do.call(max, Ms)),
#	col=cols[1]
#)
#i = 2
#for(f in files[-1]){ lines(Ms[[f]],yFills[[f]], col=cols[i]); i=i+1; }
##
#dev.new()
#plot(Ms[[1]], xFills[[1]], 'l', 
#	ylim=c(do.call(min, xFills), do.call(max, xFills)),
#	xlim=c(do.call(min, Ms), do.call(max, Ms)),
#	col=cols[1]
#)
#i = 2
#for(f in files[-1]){ lines(Ms[[f]],xFills[[f]], col=cols[i]); i=i+1; }


#files = c(
#	'bh_1200_b0_3000_5M_040_300Dat.csv', 
#	'bh_1500_b0_3000_5M_040_300Dat.csv', 
#	'bh_1800_b0_3000_5M_040_300Dat.csv', 
#	'bh_2700_b0_3000_5M_040_300Dat.csv'
#)
##files = c(
##	'bh_0900_b0_3000_M_040_300Dat.csv', 
##	'bh_1200_b0_3000_M_040_300Dat.csv', 
##	'bh_1500_b0_3000_M_040_300Dat.csv', 
##	'bh_1800_b0_3000_M_040_300Dat.csv',
##)
#Ms = list(NA, NA, NA, NA)
#names(Ms) = files
#yFills = list(NA, NA, NA, NA)
#names(yFills) = files
#xFills = list(NA, NA, NA, NA)
#names(xFills) = files

	##smaX  = SMA(dat$head.x-dat$tail.x, n=4)
	#pdf('bh_2400_b0_3000_M_040_300XBias.pdf')
	#plot(dat$M, head.x-dat$tail.x,
	#	ylab="Bias(F*/M)",
	#	xlab="Natural Mortality"
	#)
	##lines(dat$M, emaX)
	#lines(dat$M, smaX5, lwd=3)
	#lines(dat$M[1:5], fillX, lwd=3)
	##lines(dat$M, smaX, col='blue')
	###lines(dat$M, vmaX, col='green')
	##lines(dat$M, wmaX, col='green')
	##
	#dev.off()
	##
	#pdf('bh_2400_b0_3000_M_040_300YBias.pdf')
	##dev.new()
	##emaY = EMA(dat$head.y-dat$tail.y, n=4)
	##smaY  = SMA(dat$head.y-dat$tail.y, n=4)
	#plot(dat$M, head.y-dat$tail.y,
	#	ylab=expression("Bias("*B[h]/B[0]*")"),
	#	xlab="Natural Mortality"
	#)
	##lines(dat$M, emaY)
	#smaY5[1] = (head.y-dat$tail.y)[1]
	#lines(dat$M, smaY5, lwd=3)
	#lines(dat$M[1:5], fillY, lwd=3)
	##lines(dat$M, smaY, col='blue')
	###lines(dat$M, vmaX, col='green')
	##lines(dat$M, wmaY, col='green')
	#dev.off()

