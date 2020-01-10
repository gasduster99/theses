rm(list=ls())

#
library(rstan)

#
#DATA
#

#
year = 1965:1987; Y = length(year)
#specify data to be imported into stan 
D = list(#use the same R variable names here that will be used in the stan data block
	N = Y,
	cpue = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63),
	#'catch' is an internal variable to stan; no variable names are allowed to be call 'catch'
	ctch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
)

#
#FIT
#

#interface w/ stan
out = stan( file = "stanTry.stan",
	data = D,
	iter = 10^4,
	chains = 10
)
#extract samples from stan output
samples = extract(out, permuted=T)

#see posteriors
plot(out)
