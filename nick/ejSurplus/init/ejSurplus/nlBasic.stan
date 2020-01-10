data{
        int<lower=0> N;
        real<lower=0> Y[N];
	real<lower=0> X[N];
}

parameters{
	real<lower=0> a;
	real<lower=0> b;
	real<lower=0> l;
	real<lower=0> s;
}

model{
	real logMean[N];
	//vector[N] meen;
	for(i in 1:N){
		logMean[i] = log(a - b*l^X[i]);
		//meen[i] = a - b*l^X[i];
	}
	Y ~ lognormal(logMean, s);
	
}


