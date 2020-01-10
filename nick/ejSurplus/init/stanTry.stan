data{
	int<lower=0> N;
	real<lower=0> cpue[N];
	//'catch' is an internal variable to stan; no variable names are allowed to be call 'catch'
	int<lower=0> ctch[N];
}

//these get estimated
parameters{
	real<lower=0, upper=1> q;
        real<lower=0, upper=1> r;
        //real logK;
	real<lower=0> K;
        real<lower=0> s;        
}

transformed parameters{
	//real<lower=0> K;
	//K = exp(logK);
}

model{
	//build data generating model
	real B0;
	real B[N]; //no constraints here 
	real logMean[N];
	B0 = K;
	B[1] = B0 + r*B0*(1-(B0/K));
	logMean[1] = log(q) + log(B[1]);
	for(t in 2:N){
		//B[t] = B[t-1] + r*B[t-1]*(1-exp(log(B[t-1]) - logK)) - ctch[t-1];
		B[t] = B[t-1] + r*B[t-1]*(1-exp(log(B[t-1]) - log(K))) - ctch[t-1];
		logMean[t] = log(q) + log(B[t]);
	}
	//
	cpue ~ lognormal(logMean, s);
	//priors
	q ~ beta(1, 5);
	r ~ beta(3, 5);
	//logK ~ cauchy(0, 8);
	K ~ cauchy(0, 2900);
	s ~ cauchy(0, 1);
}
