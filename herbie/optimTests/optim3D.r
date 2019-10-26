rm(list=ls())

rosenbrock <- function(x){
        n <- length(x)
        out = sum(100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
        return(out)
}

rastrigin <- function(x){
        out = 10*length(x)+sum(x^2-10*cos(2*pi*x))
        return(out)
}


#h = 0.1
##h = 0.01
#xRose = seq(-5, 5, h)
#yRose = seq(-5, 5, h)
#ZRose = matrix( NaN, nrow=length(xRose), ncol=length(yRose) )
#zRast = matrix( NaN, nrow=length(xRose), ncol=length(yRose) )
#for (i in seq(1,length(xRose))){
#        for (j in seq(1,length(yRose))){
#                ZRose[i,j] = rosenbrock( c(xRose[i], yRose[j]) )
#		ZRast[i,j] = rastrigin( c(xRose[i], yRose[j]) )
#        }
#}

#top corner rosebrock
guess = c(1,0)

#optOutBFGS = optim( guess, rosenbrock, hessian=T, method='BFGS')#method="BFGS",  control=list(fnscale=-1), hessian=T )
#writeLines('\nBFGS Rosenbrock Outputs:')
#writeLines(sprintf('Initial Guess: [%2.3f, %2.3f]', guess[1], guess[2]))
#print( optOutBFGS )
#
#optOutSANN = optim( guess, rosenbrock, hessian=T, method='SANN')
#writeLines('\nSANN Rosenbrock Outputs:')
#writeLines(sprintf('Initial Guess: [%2.3f, %2.3f]', guess[1], guess[2]))
#print( optOutSANN )

optOutBFGS = optim( guess, rastrigin, hessian=T, method='BFGS')#method="BFGS",  control=list(fnscale=-1), hessian=T )
writeLines('\nBFGS Rastigin Outputs:')
writeLines(sprintf('Initial Guess: [%2.3f, %2.3f]', guess[1], guess[2]))
print( optOutBFGS )

optOutSANN = optim( guess, rastrigin, hessian=T, method='SANN')
writeLines('\nSANN Rastigin Outputs:')
writeLines(sprintf('Initial Guess: [%2.3f, %2.3f]', guess[1], guess[2]))
print( optOutSANN )
