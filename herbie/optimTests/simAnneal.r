rm(list=ls())

#a function defining the energy associate with a given state
#i.e. the function to be optimized.
E = function(ex){
	out = sin(ex) + 0.001*ex
	return(out)
}

logist = function(x,a,b){1/(1+exp(-a-b*x))}

#temperature is a stopping criteria that controls the locality of the search.
tempFunc = function(kay){
	out = const*logist(kay, a, b)+fudge
	return(out)
}

#a function for defining the acceptance probability
accept = function(guess, acpt, t){
	if (guess<acpt){
		out = 1
	}else{
		out = exp(-(guess-acpt)/(t))
	} 
	
	return(out)
}

#monitoring convergence
left = 0
right = 10*pi

#initial guess(state)
x0 = 0
#number of iterations
kMax = 600
#a constant determining how different the variance of our proposal  can get from zero (large values make the difference large, and small values make the difference small) 
const = 1.99999999
#a constant choosing the cooling rate of the process
b = 0.06 #slow
#a constant controlling the where the rapid cooling events occure
a = -b*3*kMax/4 - log(2/const - 1) #centered at 5 * 0.1^-1 = 50
fudge = 0.08

#curve(tempFunc(x), from=0, to=kMax, n=1000, main='Tempature Function')

#loading the recursive variables with the initial values
x = x0
e = E(x0)
#variables for tracking the best state and energy
xBest = x 
xBests = matrix(NaN, nrow=kMax, ncol=1)
eBest = e
eBests = matrix(NaN, nrow=kMax, ncol=1)
xS = matrix(NaN, nrow=kMax, ncol=1)
temps = matrix(NaN, nrow=kMax, ncol=1)

acc = 0
it = 0 

pdf('simPics.pdf', width=14, height=9)
layout(matrix(c(1,2,3), nrow=1, ncol=3))
curve( E(x), from=left, to=right, n=1000, main = 'f(x) = sin(x) + 0.001x', ylim=c(-1, 1.5))
for ( k in seq(1, kMax) ){
	standK = kMax-k
	temp = tempFunc(standK)
	
	#choose a new element from the domain via the proposal distribution 
	#a random walk about the domain value responsible for the current function value
	xHat = rnorm(1, x, temp) 
	#x starts in the domain, but if it ventures out of the domain try again until it stays in the domain
	while (xHat<left | xHat>right){ xHat = rnorm(1, x, temp) }	
	eHat = E(xHat)
	
	#define the accepance probability
	roe = accept(e, eHat, temp)
	#accept with probability roe
	#bernoulli p=roe
	if ( rbinom(1, 1, roe) ) {
		x = xHat
		e = eHat
		acc = acc+1
	}	

	#update our best guesses
	if ( eBest<e ){
		xBest = x
		eBest = e
	}

	#text(20, 1.25, sprintf('Temp: %f', temp))
	lines(c(x, x), c(-100, e), col='red')
	
	#get ready for the next iteration 
	temps[k] = temp
	xS[k] = x
	xBests[k] = xBest
	eBests[k] = eBest	
	it = it+1
} 

lines(c(xBest, xBest), c(-100, eBest), col='blue', lwd=2)

d=seq(1, kMax)
its = seq(1, it)

plot(its, temps, 'l', xlab='Iteration', ylab='Temperature', main='Temperature Function')

plot(its, xBests, 'l', col='green', xlab='Iteration', ylab='Value', ylim=c(0, 30), main='Convergence')
lines(its, xS, col='blue')
legend('bottomright', c('Current', 'Best'), lty=c(1,1), col=c('blue', 'green'))

dev.off()

writeLines(         'My Outputs:')
writeLines( sprintf('Guess           = %f', x0) )
writeLines( sprintf('Par             = %f', xBest) )
writeLines( sprintf('Value           = %f', eBest) )
writeLines( sprintf('Iterations      = %d', it) )
writeLines( sprintf('Acceptance Rate = %f', acc/it) )

#pdf('poly32.pdf') #main = 'f(x) = sin(x) + 0.001x'

#plot(xBests, )

#legend( 'bottomright', legend=c('f(x)', 'f\'(x)', 'f\'\'(x)'), col=c('black', 'red', 'bl
#dev.off()



#s ← s0; e ← E(s)                                  // Initial state, energy.
#sbest ← s; ebest ← e                              // Initial "best" solution
#k ← 0                                             // Energy evaluation count.
#while k < kmax and e > emax                       // While time left & not good enough:
#  T ← temperature(k/kmax)                         // Temperature calculation.
#  snew ← neighbour(s)                             // Pick some neighbour.
#  enew ← E(snew)                                  // Compute its energy.
#  if P(e, enew, T) > random() then                // Should we move to it?
#    s ← snew; e ← enew                            // Yes, change state.
#  if enew < ebest then                            // Is this a new best?
#    sbest ← snew; ebest ← enew                    // Save 'new neighbour' to 'best found'.
#  k ← k + 1                                       // One more evaluation done
#return sbest                                      // Return the best solution found.
