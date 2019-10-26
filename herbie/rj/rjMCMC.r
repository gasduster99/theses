rm(list=ls())

library(gtools)
library(MASS)
library(stats)

tSource = function(fName){
	t=proc.time()
	source(fName)
	
	return(proc.time()-t)
}

Mode = function(x) {
	ux = unique(x)
  	out = ux[which.max(tabulate(match(x, ux)))]
	
	return(out)
}

isame = function(list1, list2){
	#Give theis function two vectors (ie. list1, list2). 
	#This function will determine if the vectors are the "same".
	#I define "same", as vectors containing the same number of elements in the same order.
	#This function returns a boolean if theye are or are not the same.
	
	truth = list1==list2
	n = length(truth)
	count=0
	for (item in truth){
		if (item){count=count+1}
	}
	
	return(n==count)
}

dbk = function(km){
	#Give this function a number of components (k).
	#This function will figure out the probability of a split or merge.
	#Most of the time merge and split are equally weighted.
	#If there is {one|kMax} component the probability of mergeing is {0|1}, and splitting is thus {1|0} 
	#dk=out[1]=pr(merge) #bk=out[2]=pr(split)
	
	out = c(0.5, 0.5)
	if (km==1){out=c(0, 1)}
	if (km==kMax){out=c(1, 0)}
	if (km<1 | km>kMax){print('Run Away Dimension (k)')}
	
	return(out)
}

split = function(startList, eleStar1, eleStar2, s, last){
	#Give this function a vector(or List!), the elements that split the element that we are replacing, and the index of the split-up component.
        #This function will remove the elements at s and replace it with eleStar1 and eleStar2, in that order.
        #The merged vector is returned

	if (typeof(startList)=="list"){out=list()}else{out=c()}
	for (i in seq(1, length(startList))){	
		#the left-most or the middlest
		if ( (s==1 & i==1)|(s>1 & s<last & i==s) ){
			if (typeof(startList)=="list"){out=c(out, list(eleStar1), list(eleStar2))}else{out=c(out, eleStar1, eleStar2)}
			next
		}
		#the right-most
		if (s==last & i==(s-1)){
			if (typeof(startList)=="list"){out=c(out, list(startList[[i]]), list(eleStar1), list(eleStar2))}else{out=c(out, startList[i], eleStar1, eleStar2)}
			break
		}
		if (typeof(startList)=="list"){out=c(out, list(startList[[i]]))}else{out=c(out, startList[i])}	
	}
	
	return(out)
}			

merge = function(startList, eleStar, m1, m2, last){
	#Give this function a vector(or List! :)), a merged element, and the two adjacent indicies to merge
	#This function will remove the elements at m1 and m2 and replace then with eleStar
	#The merged vector is returned

	jLeft = min(m1, m2)
	jRight = max(m1, m2) 
	flag = 0 
	if (typeof(startList)=="list"){out=list()}else{out=c()}
	for (i in seq(1, length(startList))){	
		if (flag){flag=0;next}
		#the left-most or the middlest
		if ( (jLeft==1 & i==1)|(jLeft>1 & jRight<last & i==jLeft) ){
			if(typeof(startList)=="list"){out=c(out, list(eleStar))}else{out=c(out, eleStar)}
			flag = 1
			next
		}
		#the right-most
		if (jRight==last & i==jLeft){
			if (typeof(startList)=="list"){out=c(out, list(eleStar))}else{out=c(out, eleStar)}
			break
		}
		if (typeof(startList)=="list"){out=c(out, list(startList[[i]]))}else{out=c(out, startList[i])}	
	}
	
	return(out)
}

add = function(eleStar, startList, JStar, last){
        #the last element 
        if (sprintf('%s',JStar)=="NA"){ out=split(startList, startList[[length(startList)]], eleStar, length(startList), last)
        #first k[m]-1 elements
        }else{ out=split(startList, eleStar, startList[[JStar]], JStar, last) }

        return(out)
}

kill = function(startList, JKill, last){

	if (JKill==last){ out=merge(startList, startList[[JKill-1]], JKill, JKill-1, k[m+1])
	}else{ out=merge(startList, startList[[JKill+1]], JKill, JKill+1, k[m+1]) }

	return(out)
}

funk = function(listX, par){#, madHatter){
	#Give this function a list of vectors of data (i.e. listX) and a list of parameters (i.e. par), and an n-dimensional vector of j indices
	#The elements of par must be k[m]-dimensional

	W = par[[1]]
	MU = par[[2]]
	SIG2 = par[[3]]
	
	out = list()
	for (j in seq(1, length(listX))){
		out[[j]]=numeric(0)
		if (sprintf('%s', listX[[j]][1])=='NA'){next}
		for (i in seq(1, length(listX[[j]]))){
			x = listX[[j]][i]
			const = W/sqrt(SIG2)
			inside = -((x-MU)^2)/(2*SIG2)
			prop = const*exp(inside)
			prHat = ( prop/sum(prop) )[ j ]
			out[[j]] = c(out[[j]], prHat)
		}
	}

	return(out)	
}

funkyFunk = function(x, par){
	#Given this function a data point to evaluate at (ie. x), and a list of the k[m]-dimentional vectors W, MY, SIG2.
	#This function evaluates eq.(8) for each component at a particular observation.
	#A k[m]-dimensional vector of probabilities is returned for the given data point.   
	
	W = par[[1]]
	MU = par[[2]]
	SIG2 = par[[3]]
	const = W/sqrt(SIG2)
	inside = -((x-MU)^2)/(2*SIG2)
	prop = const*exp(inside)
	out = prop/sum(prop)	
	
	return(out)
}

allocate = function(data, par){
	#Give this function some data to be allocated, a function to determine the probabilites to determine how the allocation should occur, and the dimension to allocate into.
	#You will get allocated data(out1), the number of observations in each component(out2), and a vector of the jHats for each n in the data (i.e. out3) back.

	nn = length(data)
	
	dim = length(par[[1]])
	out1 = list()
	for (j in seq(1, dim)){out1[[j]]=numeric(0)}
	#a matrix of counters
	out2 = matrix(0, nrow=dim, ncol=1)
	#if a component is empty there is nothing to allocate
	if (nn==0){ return(list(out1, out2)) }	
	for(i in seq(1, nn)){
		pri = funkyFunk(data[i], par)	
		jHat = sample(seq(1, dim), 1, replace=T, prob=pri)
		out1[[jHat]] = c(out1[[jHat]], data[i])
		out2[jHat] = out2[jHat]+1
	}

	return( list(out1, out2) )
}

nEmpty = function(weye){
	dim = length(weye)
        out = matrix(0, nrow=dim, ncol=1)
        for(j in seq(1, length(weye))){ if(length(weye[[j]])==0){out[j]=1} }

        return(out)
}

organize = function(mixedList){
	#a function to figure out all of the samples with j components then make a list where the number of components are in the list indices
	#the elements of the list are matricies with j columns and n_j rows
        
	out = list()
        for (j in seq(1, max(k))){out[[j]]=numeric(0)}
        for(ind in seq(1, length(mixedList))){
                jHat = length(mixedList[[ind]])
                out[[jHat]] = rbind(out[[jHat]], mixedList[[ind]])
        }   
        
	return(out)
}

likelyhood = function(data, par){
	#Give this function 2 lists, data is a list containing the data split up by component.
	#par is a list containnig the parameters needed for a mixture of normals.
	#each element of the par list is a k[m]-dimensional vector for each respective parameter. 
	
	w = par[[1]]
	mu = par[[2]]
	s2 = par[[3]]
	
	out = 1
	#given z, no w[j]
	for (j in seq(1, length(data))){ out=out*prod(w[j]*dnorm(data[[j]], mu[j], sqrt(s2[j]))) }
	
	return(log(out))
}

prior = function(uu, mu, s2, beta, k){
	#puu = ddirichlet(uu, rep(delta, k)) 
	#pmu = prod(dnorm(mu, squig, 1/sqrt(kappa)))
	#overPs2 = prod(dgamma(s2, alpha, beta))
	#ps2 = 1/overPs2
	#pz = prod(uu) 
	#pbeta = dgamma(beta, shape=g, rate=h)
	#pk = dpois(k, lambda)	

	#out = puu * pmu * ps2 * pk * pbeta# * pz# * pbeta

	lpuu = log(ddirichlet(uu, rep(delta, k))) 
	lpmu = sum(dnorm(mu, squig, 1/sqrt(kappa), log=T))
	loverPs2 = sum(dgamma(s2, alpha, beta, log=T))
	lps2 = -loverPs2
	lpz = sum(log(uu)) 
	lpbeta = dgamma(beta, shape=g, rate=h, log=T)
	lpk = dpois(k, lambda, log=T)

	lout = (lpuu + lpmu + lps2 + lpz + lpbeta + lpk)

	#writeLines(sprintf('lpuu:%f\nlpmu:%f\nloverPs2:%f\nlps2:%f\nlpz:%f\nlpbeta:%f\nlpk:%f\nlout:%f\n', lpuu, lpmu, loverPs2, lps2, lpz, lpbeta, lpk, lout))
	#writeLines(sprintf('puu:%f\npmu:%f\noverPs2:%f\nps2:%f\npz:%f\npbeta:%f\npk:%f\nout:%f\n', puu, pmu, overPs2, ps2, pz, pbeta, pk, out))

	return(lout)
}


#grades for ams206 winter 2012
classAVG = c(0.5692307692,0.5346153846,0.8615384615,0.8923076923,0.6192307692,0.7384615385,0.8384615385,0.6884615385,0.6730769231,0.7192307692, 0.63846153, 0.85,0.8269230769,0.9346153846,0.9230769231,0.6076923077,0.9,0.4961538462,0.9, 0.76923077,0.7615384615,0.8807692308,0.6423076923,0.6961538462,0.5211538462,    0.5153846154,0.9423076923,0.55,0.55)
n = length(classAVG)

#RANGE
lower = min(classAVG)
upper = max(classAVG)
R = upper-lower

#PARAMETERS
squig = lower+(R/2)
kappa = 1/(R^2) #a small multiple
ks = kappa*squig
alpha = 1.5
g = 0.5
h = 1/(R^2) #a small multiple
lambda = 1
delta = 1

#CONVERGENCE
kMin = 1
kMax = 5#length(classAVG)#maybe just 10
maxIt = 10000#100000#

#DATA STRUCTURES
beta = matrix(NaN, nrow=maxIt+1, ncol=1)
k = matrix(NaN, nrow=maxIt+1, ncol=1)
muRJ = list()
sig2RJ = list()
wRJ = list()
likelyHood = matrix(NaN, nrow=maxIt+1, ncol=1)
propPost = matrix(NaN, nrow=maxIt+1, ncol=1)

#INTIALIZING
z = kMin
beta[1] = 1
k[1] = kMin
muRJ = c(muRJ, lower + (R/2))
sig2RJ = c(sig2RJ, 1)
wRJ = c(wRJ, 1/kMin)
likelyHood[1] = likelyhood(list(classAVG), list(wRJ[[1]], muRJ[[1]], sig2RJ[[1]]))
propPost[1] = likelyHood[1] * prior(wRJ[[1]], muRJ[[1]], sig2RJ[[1]], beta[1], k[1])

mCount = 0
maCount = 0

sCount = 0
saCount = 0

dCount = 0
daCount = 0

bCount = 0
baCount = 0

#ITERATING
#Y is a list of the observation that belong to each component, that recursively updates
#I need to make a list so that I can have a row for each component             
Y = list(classAVG)
#ns is the number of elements in each of the rows of Y
ns = matrix(length(Y[[k[1]]]), nrow=k[1], ncol=1)
for (m in seq(1, maxIt)){	
	#WEIGHTS
        wPVect = delta + ns  
        wDraw = c(rdirichlet(1, wPVect))
        wRJ = c(wRJ, list(wDraw))
	
	#MEANS	
	muP1 = c()
	for (j in seq(1,k[m])){ 
		muP1[j] = ( (sum(Y[[j]])/sig2RJ[[m]][j]) + ks )/( (ns[j]/sig2RJ[[m]][j]) + kappa )
	}
	#vectorized to make and independent normal covariance matrix
	muP2 = sqrt( (ns/sig2RJ[[m]] + kappa)^(-1) )
	muDraw = rnorm(k[m], muP1, muP2)
	#keeping it ordered
	if (isame(muDraw, sort(muDraw))){muRJ=c(muRJ, list(muDraw))}else{muRJ=c(muRJ, list(muRJ[[m]]))}

	#VARIANCES
	sig2 = c()
	for (j in seq(1,k[m])){
		a = alpha + ns[j]/2 
		b = beta[m] + sum( (Y[[j]] - muRJ[[m+1]][j])^2 )/2
		overSig2 = rgamma(1, shape=a, rate=b)
		sig2[j] = 1/overSig2
	}
	sig2RJ = c(sig2RJ, list(sig2))

	#ALLOCATION
	par = list(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]])
	yn = allocate(classAVG, par)	
	Y = yn[[1]]
	ns = yn[[2]]

	#BETA
	a = g + k[m]*alpha
	b = h + sum( 1/sig2RJ[[m+1]] )
	beta[m+1] = rgamma(1, shape=a, rate=b)
	
	k[m+1] = k[m]

	#MERGE/SPLIT
	#dk=db[1]=pr(merge) #bk=db[2]=pr(split) 
	db = dbk(k[m])
	 	
	#MERGE
	if (rbinom(1, 1, db[1])){
		#merge points
		jm1 = sample(seq(1, k[m]), 1)	
		#gotta consider the case where you merge the first or last component of the current model, else you merge some components in the middle
		if(jm1 == 1){jm2=2}else if(jm1 == k[m]){jm2=k[m]-1}else{jm2=sample(c(jm1-1, jm1+1), 1)}
	
		#WEIGHTS(merge)
		wStar = wRJ[[m+1]][jm1] + wRJ[[m+1]][jm2] 

		#MEANS(merge)
		one = wRJ[[m+1]][jm1]*muRJ[[m+1]][jm1]
		two = wRJ[[m+1]][jm2]*muRJ[[m+1]][jm2]		
		muStar = (one + two)/wStar 

		#VARIANCES(merge)
		one = wRJ[[m+1]][jm1]*(muRJ[[m+1]][jm1]^2 + sig2RJ[[m+1]][jm1]^2)
		two = wRJ[[m+1]][jm2]*(muRJ[[m+1]][jm2]^2 + sig2RJ[[m+1]][jm2]^2)
		sig2Star = (one + two)/wStar - muStar^2

		#REALLOCATE(merge)		
		YStar = c(Y[[jm1]], Y[[jm2]])
	
		wMerge = merge(wRJ[[m+1]], wStar, jm1, jm2, k[m])
		muMerge =  merge(muRJ[[m+1]], muStar, jm1, jm2, k[m])
                sig2Merge = merge(sig2RJ[[m+1]], sig2Star, jm1, jm2, k[m])
		YMerge = merge(Y, YStar, jm1, jm2, k[m])
		nsMerge = merge(ns, ns[jm1]+ns[jm2], jm1, jm2, k[m])

		par = list(wMerge, muMerge, sig2Merge)		
		prAllMerge = funk(YMerge, par)
		prAll = prod( unlist(prAllMerge) )		
		
		#likeRatio = 1
		logLikeRatio = 0     
                #likeNow                                                    log(wRJ[[m+1]][j])+                     
                for (j in seq(1, length(Y))){logLikeRatio=logLikeRatio-sum( log(wRJ[[m+1]][j])+dnorm(Y[[j]], muRJ[[m+1]][j], sqrt(sig2RJ[[m+1]][j]), log=T) )}
                #likeMerge                                                       log(wMerge[j])+
                for (j in seq(1, length(YMerge))){logLikeRatio=logLikeRatio+sum( log(wMerge[j])+dnorm(YMerge[[j]], muMerge[j], sqrt(sig2Merge[j]), log=T) )}
		
		l1 = ns[jm1]
		l2 = ns[jm2]
		
		u1 = rbeta(1, 2, 2)
                u2 = rbeta(1, 2, 2)
		u3 = rbeta(1, 1, 1)
		
		#PR(ACCEPT)
		logq1 = logLikeRatio + dpois(k[m]-1, lambda, log=T) - dpois(k[m], lambda, log=T)
		logLineOne = logq1 + log(1/k[m]) + log(wStar^(delta-1+l1+l2)) - log(wRJ[[m+1]][jm1]^(delta-1+l1)*wRJ[[m+1]][jm2]^(delta-1+l2)*beta(delta, k[m]*delta))
		logLineTwo = log(sqrt(kappa/(2*pi))) + ( (-kappa/2)*((muStar-squig)^2 - (muRJ[[m+1]][jm1]-squig)^2 - (muRJ[[m+1]][jm2]-squig)^2) )
		logLastBit = (-beta[m+1]*(1/sig2Star - 1/sig2RJ[[m+1]][jm1] - 1/sig2RJ[[m+1]][jm2])) 
		logLineThree = log(beta[m+1]^alpha) - lgamma(alpha) + (-alpha-1)*log(sig2Star/(sig2RJ[[m+1]][jm1]*sig2RJ[[m+1]][jm2])) + logLastBit
		logq2 = log(dbk(k[m]-1)[1]) - log(db[2]) -log(prAll)
		logLineFour = logq2 + log(dbeta(u1, 2, 2)*dbeta(u2, 2, 2)*dbeta(u3, 1, 1)) 		
		logLineFive = log(u2*(1-u2^2)*u3*(1-u3)*sig2Star) - log(wStar * abs(muRJ[[m+1]][jm1]-muRJ[[m+1]][jm2]) * sig2RJ[[m+1]][jm1] * sig2RJ[[m+1]][jm2])
		A = exp(logLineOne + logLineTwo + logLineThree + logLineFour + logLineFive)
		
		mCount = mCount + 1 
		
		#ROLL DICE
		roe = min(1, A)
		if (rbinom(1,1,roe)){
			#ACCEPT
			#print('MERGE')	
			wRJ[[m+1]] = wMerge
			muRJ[[m+1]] = muMerge
			sig2RJ[[m+1]] = sig2Merge	
			Y = YMerge
			ns = nsMerge
			k[m+1] = k[m]-1
			maCount = maCount + 1
		}	
	#SPLIT	
	}else{
		#split points
		js = sample(seq(1, k[m]), 1)	
		u1 = rbeta(1, 2, 2)
		u2 = rbeta(1, 2, 2)
		u3 = rbeta(1, 1, 1)

		#WEIGHTS(split)
		wStar1 = wRJ[[m+1]][js]*u1 
		wStar2 = wRJ[[m+1]][js]*(1-u1)
		wSTAR = c(wStar1, wStar2)
		
		#MEANS(split)
		muStar1 = muRJ[[m+1]][js] - u2*sqrt(sig2RJ[[m+1]][js]*(wStar2/wStar1))
		muStar2 = muRJ[[m+1]][js] + u2*sqrt(sig2RJ[[m+1]][js]*(wStar1/wStar2))
		muSTAR = c(muStar1, muStar2)
		#ADJACENT?
		minny = which(muSTAR==min(muSTAR))
		maxxy = which(muSTAR==max(muSTAR))
		muSplit = split(muRJ[[m+1]], muSTAR[minny], muSTAR[maxxy], js, k[m])
		#the new mu's aren't adjacent, so you automatically reject
		if (isame(muSplit, sort(muSplit))==F){
			likelyHood[m+1] = likelyhood(Y, list(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]]))
			propPost[m+1] = likelyHood[m+1] * prior(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]], beta[m+1], k[m+1])
			next
		}
		
		#VARIANCES(split)
		sig2Star1 = u3*(1-u2^2) * sig2RJ[[m+1]][js] * (wRJ[[m+1]][js]/wStar1)
		sig2Star2 = (1-u3)*(1-u2^2) * sig2RJ[[m+1]][js] * (wRJ[[m+1]][js]/wStar2)
		sig2STAR = c(sig2Star1, sig2Star2)
		
		#REALLOCATE
		par = list(wSTAR, muSTAR, sig2STAR)
		yn = allocate(Y[[js]], par)
		YSTAR = yn[[1]]
		nsSTAR = yn[[2]]
		
		l1 = nsSTAR[1]
		l2 = nsSTAR[2]
	
		#muSplit above, MEANS(split) section
		wSplit = split(wRJ[[m+1]], wSTAR[minny], wSTAR[maxxy], js, k[m]) 
		sig2Split = split(sig2RJ[[m+1]], sig2STAR[minny], sig2STAR[maxxy], js, k[m]) 
		YSplit = split(Y, YSTAR[[minny]], YSTAR[[maxxy]], js, k[m])#[[1]] [[2]]
		nsSplit = split(ns, nsSTAR[minny], nsSTAR[maxxy], js, k[m])#[[1]] [[2]]

		par = list(wSplit, muSplit, sig2Split)
		prAllSplit = funk(YSplit, par)
		prAll = prod(unlist(prAllSplit))

		likeRatio = 1		
		#likeNow                                        *      wRJ[[m+1]][j]* 
		for (j in seq(1, length(Y))){likeRatio=likeRatio/prod( wRJ[[m+1]][j]*dnorm(Y[[j]], muRJ[[m+1]][j], sqrt(sig2RJ[[m+1]][j])) )}
		#likeSplit                                           /      wSplit[j]* 
		for (j in seq(1, length(YSplit))){likeRatio=likeRatio*prod( wSplit[j]*dnorm(YSplit[[j]], muSplit[j], sqrt(sig2Split[j])) )}

		#Pr(ACCEPT)
		q1 = likeRatio * (dpois(k[m]+1, lambda)/dpois(k[m], lambda)) * (k[m]+1) 
		lineOne = q1 * (k[m]+1) * (wStar1^(delta-1+l1)*wStar2^(delta-1+l2))/(wRJ[[m+1]][js]^(delta-1+l1+l2)*beta(delta, k[m]*delta))
		lineTwo = sqrt(kappa/(2*pi)) * exp( (-kappa/2)*((muStar1-squig)^2 + (muStar2-squig)^2 - (muRJ[[m+1]][js]-squig)^2) )
		lineThree = (beta[m+1]^alpha)/(gamma(alpha)) * (sig2Star1*sig2Star2/sig2RJ[[m+1]][js])^(-alpha-1) * exp(-beta[m+1]*(1/sig2Star1 + 1/sig2Star2 - 1/sig2RJ[[m+1]][js]))
		q2 = dbk(k[m]+1)[1]/(db[2] * prAll)
		lineFour = q2 * 1/(dbeta(u1, 2, 2)*dbeta(u2, 2, 2)*dbeta(u3, 1, 1))
		lineFive = (wRJ[[m+1]][js]*abs(muStar1-muStar2)*sig2Star1*sig2Star2)/(u2*(1-u2^2)*u3*(1-u3)*sig2RJ[[m+1]][js])
		A = lineOne * lineTwo * lineThree * lineFour * lineFive
		#writeLines(sprintf('q1:%f\nlineOne:%f\nlineTwo:%f\nlineThree:%f\nq2:%f\nlineFour:%f\nlineFive:%f\nA:%f', q1, lineOne, lineTwo, lineThree, q2, lineFour, lineFive, A))

		sCount = sCount + 1

		#ROLL DICE
		roe = min(1, A)
		if (rbinom(1,1,roe)){
			#ACCEPT
			#print('SPLIT')	
			wRJ[[m+1]] = wSplit        
			muRJ[[m+1]] = muSplit      
			sig2RJ[[m+1]] = sig2Split  
			Y = YSplit
			ns = nsSplit
			k[m+1] = k[m]+1	
			saCount = saCount +1
		}
	}
	

	#BIRTH/DEATH
        #dk=db[1]=pr(death) #bk=db[2]=pr(birth) 
        db = dbk(k[m+1])

        #DEATH
        if (rbinom(1, 1, db[1])){
		fillState = nEmpty(Y)
		k0 = sum(fillState)

		#only kill if there is an empty component	
		if (k0==0){
			likelyHood[m+1] = likelyhood(Y, list(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]]))
			propPost[m+1] = likelyHood[m+1] * prior(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]], beta[m+1], k[m+1])
			next
		}

		poss = which(fillState==1)
		jKill = sample(poss, 1)

		wStar = wRJ[[m+1]][jKill] 

		wDeath = kill(wRJ[[m+1]], jKill, k[m+1])
		muDeath = kill(muRJ[[m+1]], jKill, k[m+1])
		sig2Death = kill(sig2RJ[[m+1]], jKill, k[m+1])
		YDeath = kill(Y, jKill, k[m+1])
		nsDeath = kill(ns, jKill, k[m+1])
		
		one = (dpois(k[m+1]-1, lambda)/dpois(k[m+1], lambda)) 
		two = (wStar^(delta-1) * (1-wStar)^(n + k[m+1]*delta - k[m+1])) / (beta(k[m+1]*delta, delta) * k[m+1])
		three = (dbk(k[m+1]-1)[1] * (1-wStar)^(k[m+1])) / (db[2] * (k0+1) * dbeta(wStar, 1, k[m+1]))
		A = one * two * three	

		dCount = dCount+1

                #ROLL DICE
                roe = min(1, A)
                if (rbinom(1,1,roe)){
			#ACCEPT 
                        #print('DEATH') 
                        wRJ[[m+1]] = wDeath
                        muRJ[[m+1]] = muDeath
                        sig2RJ[[m+1]] = sig2Death
                        Y = YDeath
                        ns = nsDeath
                        k[m+1] = k[m+1]-1
                        daCount = daCount+1
                }	
	
	#BIRTH
	}else{
		#WEIGHTS(birth)
		wStar = rbeta(1, 1, k[m+1])
		#wBirth + wStar should now sum to 1
                wB = wRJ[[m+1]]*(1-wStar)
		
		#MEANS(birth)
		muStar = rnorm(1, squig, 1/sqrt(kappa))
		#if muStar is the biggest this NA's, but I use the NA in add as a flag to mean the last index
		jStar = which(muRJ[[m+1]]>muStar)[1]	
		
		#VARIANCES(birth)
		overSig2Star = rgamma(1, shape=alpha, rate=beta[m+1])
		sig2Star = 1/overSig2Star 

		wBirth = add(wStar, wB, jStar, k[m+1])
		muBirth = add(muStar, muRJ[[m+1]], jStar, k[m+1])	
		sig2Birth = add(sig2Star, sig2RJ[[m+1]], jStar, k[m+1])
		YBirth = add(numeric(0), Y, jStar, k[m+1])
		nsBirth = add(0, ns, jStar, k[m+1])

		k0 = sum(nEmpty(Y))	

		one = (dpois(k[m+1]+1, lambda)/dpois(k[m+1], lambda))
		two = (wStar^(delta-1) * (1-wStar)^(n + k[m+1]*delta - k[m+1]) * (k[m+1]+1))/beta(k[m+1]*delta, delta)
		three = (dbk(k[m+1]+1)[1] * (1-wStar)^(k[m+1])) / (db[2] * (k0+1) * dbeta(wStar, 1, k[m+1])) 
		A = one * two * three
		#writeLines(sprintf('one:%f\ntwo:%f\nthree:%f\nA:%f', one, two, three, A))

		bCount = bCount+1

                #ROLL DICE
                roe = min(1, A)
                if (rbinom(1,1,roe)){
			#ACCEPT
                    	#print('BIRTH')  
                        wRJ[[m+1]] = wBirth
                        muRJ[[m+1]] = muBirth
                        sig2RJ[[m+1]] = sig2Birth
                        Y = YBirth
                        ns = nsBirth
                        k[m+1] = k[m+1]+1
			baCount = baCount+1
                }

	}
	
	#evaluate likelyhood and propPosterior
	likelyHood[m+1] = likelyhood(Y, list(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]]))
	propPost[m+1] = likelyHood[m+1] * prior(wRJ[[m+1]], muRJ[[m+1]], sig2RJ[[m+1]], beta[m+1], k[m+1])
}


ma = maCount/mCount
sa = saCount/sCount
da = daCount/dCount
ba = baCount/bCount

wO = organize(wRJ)
muO = organize(muRJ)
sig2O = organize(sig2RJ)

system('cvlc ~/Music/samples/DING2.mp3 &')
save.image(sprintf('rjMCMC%d.RData', maxIt))








