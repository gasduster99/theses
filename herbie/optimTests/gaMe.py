import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import numpy.random as nrd
from time import sleep

'''
PARAMETERIZATION
'''
def f(x):
        out = np.sin(x) + 0.001*x
        return(out)

#must be an even number
precision = 2
popSize = 50
maxIt = 50
left = 0
right = 4*np.pi
R100 = round(right, 2)*100
#the probability that there will be a cross over site at any given bp
prCross = 0.2


'''
INITIALIZING THE FIRST POPULATION
'''
bestEst = []
meanEst = []

domain = np.arange(left, right+0.01, 0.01)
y = []
for d in domain:
	y.append( f(d) )	

fig = plt.figure( figsize=(14, 9), dpi=100, facecolor='w' )
ax = plt.subplot(5,10,1)
ax.set_xticklabels([])
plt.plot(domain, y, color='k')
plt.ylabel('f(x)')
plt.title('1')

fitness = []
#a decimal list of the population to 2 decimal places 
popDec = []
#a decimal list of the population to 2 decimal places and multiplied by 100
popDec100 = []
#a binary list of the population to 2 decimal places and multiplied by 100
popBin100 = []
for i in range(popSize):
	mem = round(rnd.choice(domain), 2)
	mem100 = int(mem*100)
	popDec.append( mem )
	popDec100.append( mem100 )
	#convert to binary then remove the leading 'ob'
	popBin100.append( bin(mem100)[2:] )

	fit = f(mem)
	fitness.append(fit)
maxEst = popDec[fitness.index(max(fitness))]
bestEst.append( maxEst/maxEst )
meanEst.append( np.mean(popDec)/maxEst )
plt.scatter(popDec, fitness, color='r')#nrd.rand(3,1))


'''
MATE AND ITERATE 
'''
for it in range(maxIt-1):
	M = max(fitness)
	m = min(fitness)
	
	I = fitness.index(M)
	i = fitness.index(m)
	
	#recplacing the least fit individual with the most fit individual
	popDec[i] = popDec[I]
	popDec100[i] = popDec100[I]
	popBin100[i] = popBin100[I]
	
	# a random vector of indecies of the first half of the population
	rndInx1 = rnd.sample( range(popSize/2), popSize/2 )
	#the second half
	rndInx2 = rnd.sample( range(popSize/2, popSize), popSize/2 )

	newBin100 = []
	newDec100 = []
	newDec = []
	newFitness = []
	#now cross-over 
	for m in range(popSize/2):
		#find the mates
		mate1 = popBin100[rndInx1[m]]
		mate2 = popBin100[rndInx2[m]]
		mates = [mate1, mate2]
		matesLens = [len(mate1), len(mate2)]
		rng = min( matesLens )
		rngI = matesLens.index(rng)
		refMate = list(mates[rngI]) 
		otherMate = list(mates[1-rngI])
	
		#a random number of cross-overs at the following indicies
		crossLoc = []
		for j in range(rng): 
			if nrd.binomial(1, prCross): crossLoc.append(j)
	
		for loc in crossLoc:
			newRef = refMate[:loc] + otherMate[loc:]
			newOther = otherMate[:loc] + refMate[loc:]
			refMate = newRef
			otherMate = newOther
		
		bin100Ref = ''.join(refMate)
		bin100Other =  ''.join(otherMate)
		dec100Ref = int(bin100Ref, 2)
		dec100Other = int(bin100Other, 2)
		#if you've wandered out of the domain then mutate until you get all the points in the domain
		while dec100Ref>R100 or dec100Other>R100:
			if dec100Ref>R100:
				#pick a random bit then flip it
				mutLoc = rnd.randrange(len(refMate))
				refMate[mutLoc] = str(1-int(refMate[mutLoc]))
				bin100Ref = ''.join(refMate)
				dec100Ref = int(bin100Ref, 2)

			else:
				#pick a random bit then flip it
                                mutLoc = rnd.randrange(len(otherMate))
				otherMate[mutLoc] = str(1-int(otherMate[mutLoc]))
				bin100Other = ''.join(otherMate)
				dec100Other = int(bin100Other, 2)

		fitRef = f(dec100Ref/100.)
		fitOther = f(dec100Other/100.)
		newBin100.append( bin100Ref )
		newBin100.append( bin100Other )
		newDec100.append( dec100Ref )
		newDec100.append( dec100Other )
		newDec.append( dec100Ref/100. )
		newDec.append( dec100Other/100. )
		newFitness.append( fitRef )
		newFitness.append( fitOther )

	popDec = newDec
	popDec100 = newDec100
	popBin100 = newBin100
	fitness = newFitness
	
	ax = plt.subplot(5, 10, it+2)
	ax.set_xticklabels([])
	if (it+2)%10!=1:ax.set_yticklabels([])
	plt.plot(domain, y, color='k')
	plt.scatter(popDec, fitness, color='r')#nrd.rand(3,1))
	if it>38: plt.xlabel('x')
	if (it+2)%10==1: plt.ylabel('f(x)')
	plt.title('%d'%(it+2))

	maxEst = popDec[fitness.index(max(fitness))]
	bestEst.append( maxEst/maxEst )
	meanEst.append( np.mean(popDec)/maxEst )

plt.tight_layout()
fig.savefig('gaSearchingMe.pdf')

maxEst = max(fitness)
argMaxEst = popDec[fitness.index(max(fitness))] 

print('      My Genetic Algorithm')
print('________________________________\n')
print('Population Size       = %d'% popSize)
print('Number of Generations = %d'% maxIt)
print('Crossover Probability = %1.2f'% prCross)
print('Search Domain: [%2.2f, %2.2f]'% (left, right) )
print('Precision             = %d\n'%precision)
print('            Results  ')
print('________________________________\n')
print('Fitness               = %2.6f'%maxEst)
print('argmax                = %2.2f\n'%argMaxEst)

fig = plt.figure( figsize=(8, 7), dpi=100, facecolor='w' )
p1, = plt.plot(range(1,maxIt+1), bestEst, color='g')
p2, = plt.plot(range(1,maxIt+1), meanEst, color='b')
plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.legend([p1, p2], ['Best', 'Mean'], loc='lower right')
plt.ylim(0, 1.1)
fig.savefig('gaConvergeMe.pdf')



























	
'''
	print(''.join(refMate), ''.join(otherMate))

	#this only produces a single cross-over site 
	for j in range(rng):
		#rbern for the probability of a crossover event
		if nrd.binomial(1, prCross):
			print('cross-over @ %d'%j)
			recomb1[j:] = otherMate[j:]
			recomb2[j:] = refMate[j:] 			

	#print(''.join(refMate), ''.join(otherMate))
	print('\n')
	newPopBin100.append( ''.join(recomb1) )
	newPopBin100.append( ''.join(recomb2) )

popBin100 = newPopBin100
'''



'''
for i in range(maxIt):
	fitness = []
	for m in range(popSize):
		fit = f(popDec[m])
		fitness.append(fit)
		
'''		
















'''
print(popDec)
print('')
print(popDec100)
print('')
print(popBin100)
'''



'''
fEval = []
	for xi in domain:
		foo = f(xi)
		fEval.append(foo)

	plt.plot(domain, fEval)
plt.show()
'''



'''
bin(10)
int(bin(10), 2)
'''
