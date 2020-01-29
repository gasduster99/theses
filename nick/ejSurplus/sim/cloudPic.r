rm(list=ls())

#
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
#DATA
#

#
n = 10000
#
l0 = cbind(rnorm(n, 10, 100), rnorm(n, 10, 100))
l1 = cbind(rnorm(n, 20, 10), rnorm(10000, 20, 10))
l2 = cbind(rnorm(n, -20, 10), rnorm(10000, 20, 10))
l3 = cbind(rnorm(n, 20, 10), rnorm(10000, -20, 10))
l4 = cbind(rnorm(n, -20, 10), rnorm(10000, -20, 10))

#
#PLOT
#

#
dev.new()
cols = brewer.pal(10, 'Spectral')
smoothScatter(l0, 
	xlim = c(min(l0[,1]), max(l0[,1])), xaxs="i",
	ylim = c(min(l0[,2]), max(l0[,2])), yaxs="i",
	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[1], 0.9)), alpha=T),
	nrpoints = 0
)
smoothScatter(l1, 
	add=T,
	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[2], 0.9)), alpha=T),
	nrpoints = 0
)
smoothScatter(l2, 
	add=T,
	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[3], 0.9)), alpha=T),
	nrpoints = 0
)
smoothScatter(l3, 
	add=T,
	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[4], 0.9)), alpha=T),
	nrpoints = 0
)
smoothScatter(l4, 
	add=T,
	colramp = colorRampPalette(c(addalpha('white', 0), addalpha(cols[5], 0.9)), alpha=T),
	nrpoints = 0
)

