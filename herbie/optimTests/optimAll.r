rm(list=ls())

f = function(x){
	out = x^3 + x^2#sin(x) + 0.001*x
        return(out)
}

d_dx = function(fun, x0, h){
	#fun is a function to take the numerical derivative 
	#x0 is the stop on the domain where you should evaulate you derivative at
	#h is the step size to use to take the numerical derivative (should be small)
	out = ( fun(x0+h)-fun(x0) )/h
	return(out)
}

d2_dx2 = function(fun, x0, h){
        #fun is a function to take the numerical derivative 
        #x0 is the stop on the domain where you should evaulate you derivative at
        #h is the step size to use to take the numerical derivative (should be small)
	fPrime = function(x0){
		out1 = d_dx(fun, x0, h)
		return(out1)
	}

	out2 = ( fPrime(x0+h)-fPrime(x0) )/h
	return(out2)
}


#smaller number start to break the algorithm due to underflow/overflow problems?
step = 10^(-6)
left = -2
right = 1

#pdf('poly32.pdf') #main = 'f(x) = sin(x) + 0.001x'
#curve( f(x), from=left, to=right, n=1000, main = 'f(x) = x^3 + x^2', ylim=c(-1, 1.5))
#curve( d_dx(f, x, step), from=left, to=right, n=1000, add=T, col='red' )
#curve( d2_dx2(f, x, step), from=left, to=right, n=1000, add=T, col='blue' )
#legend( 'bottomright', legend=c('f(x)', 'f\'(x)', 'f\'\'(x)'), col=c('black', 'red', 'blue'), lty=c(1,1,1) )
#dev.off()

maxI = 50
its = 1:maxI 
initGuess = -2
guess = initGuess
guessesNeg2 = matrix(NaN, nrow=maxI, ncol=1)
guessesNeg2[1,1] = guess
for (i in its[2:length(its)]){
	frac = d_dx(f, guess, step) / d2_dx2(f, guess, step)
	guess = guess - frac
	guessesNeg2[i, 1] = guess
}

initGuess = 0
guess = initGuess
guesses0 = matrix(NaN, nrow=maxI, ncol=1)
guesses0[1,1] = guess
for (i in its[2:length(its)]){
	frac = d_dx(f, guess, step) / d2_dx2(f, guess, step)
	guess = guess - frac
	guesses0[i, 1] = guess
}

initGuess = 2
guess = initGuess
guesses2 = matrix(NaN, nrow=maxI, ncol=1)
guesses2[1,1] = guess
for (i in its[2:length(its)]){
	frac = d_dx(f, guess, step) / d2_dx2(f, guess, step)
	guess = guess - frac
	guesses2[i, 1] = guess
}

pdf('gradConvergePoly32.pdf')
plot(its, guessesNeg2[its, 1], 'l', xlab='Iterations', ylab='Best Guess', ylim=c(-2,2))
lines(its, guesses0[its, 1], col='red', lty=2)
lines(its, guesses2[its, 1], col='blue')
legend('bottomright', title='Start at:', legend=c('-2', '0', '2'), col=c('black', 'red', 'blue'), lty=c(1,2,1))
dev.off()

#writeLines(         'My Outputs:')
#writeLines( sprintf('Guess      = %f', initGuess) )
#writeLines( sprintf('Par        = %f', guess) )
#writeLines( sprintf('Value      = %f', f(guess)) )
#writeLines( sprintf('Iterations = %d', maxI) )
#writeLines( sprintf('Hessian    = %f', d2_dx2(f, guess, step)) )
#
##optim minimizes, but fnscale negative turns it into a maximum
#optOut = optim( guess, f, method="BFGS",  control=list(fnscale=-1), hessian=T )
#writeLines('\nOptim Outputs:')
#print( optOut )
