rm(list=ls())

library(rstan)

#data for stan
D = list(
	N = 27,
        X = c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10, 12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5),
	Y = c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.7, 2.72, 2.57)
)
#
threads = 8
#interface w/ stan
out = stan( file = "nlBasic.stan",
        data = D,
        iter = 10^4,
	#control=list(
	#	adapt_delta=0.7
	#),
        chains = threads*10,
	cores = threads
)

