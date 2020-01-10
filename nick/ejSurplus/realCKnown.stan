data{
	int<lower=0> N;
	real<lower=0> cpue[N];
	int<lower=0> ctch[N];
	int<lower=0> k;
}

//these get estimated
parameters{
	real<lower=0, upper=1> q;
        real<lower=0, upper=1> sigma;
	real<lower=0> lambda; //? bounds
	real<lower=0> kappa;  //? bounds
	real<lower=0> alpha;  //? bounds
	real<lower=0> beta;   //? bounds
	real<lower=0> gamma;  //? bounds
	real<lower=0> s;
}

transformed parameters{
}

model{
	//data model	
	real logR[N];
	real B[N];
	real W[N]; 
	real logI[N];
	real wk;
	int tk;
	int C0;
	real R0;
	real B0;
	real W0;	
	//
	R0 = 1;    //?value
	B0 = 2900; //?value
	W0 = 3000; //?value
	C0 = 0;    //?value
	//
	wk = 0;    //? which value
	for(a in 1:k){ 
		wk = lambda + kappa*wk;
	}
	//
	B[1] = R0 + sigma * (lambda+kappa*W0) * (B0-C0) / W0;
	logR[1] = log(alpha) + log(B0-C0) + log(1-beta*gamma*(B0-C0))/gamma;
	W[1] = B[1]/( exp(logR[1])/wk + sigma*(B0-C0)/W0 );
	logI[1] = log(q) + log(B[1]-ctch[1]/2.);	
	//
	for(t in 2:N){	
		//
		if(t>k){
			//
			tk = t-k;
			logR[t] = log(alpha) + log(B[tk]-ctch[tk]) + log(1-beta*gamma*(B[tk]-ctch[tk]))/gamma;
		}else{
			logR[t] = log(R0);
		}
		//
		B[t] = exp(logR[t]) + sigma * (lambda+kappa*W[t-1]) * (B[t-1]-ctch[t-1]) / W[t-1];
		W[t] = B[t]/( exp(logR[t])/wk + sigma*(B[t-1]-ctch[t-1])/W[t-1] );
		logI[t] = log(q) + log(B[t]-ctch[t]/2.);
	}
	cpue ~ lognormal(logI, s);
	
	//parameter model
        q ~ beta(1, 5);
	sigma ~ beta(1, 1); //? values
	lambda ~ cauchy(0, 1);
	kappa ~ cauchy(0, 1);
	alpha ~ cauchy(0, 1);
	beta ~ cauchy(0, 1);
	gamma ~ cauchy(0, 1);
        s ~ cauchy(0, 1);
	

	//real R[N];
	//R[1] = alpha * (B[t-k]-ctch[t-k]) * (1-beta*gamma*(B[t-k]-ctch[t-k]))^(1/gamma);
	//
	//logMean[1] = log(q) + log(B[1]);
	//for(t in 2:N){
	//	//B[t] = B[t-1] + r*B[t-1]*(1-exp(log(B[t-1]) - logK)) - ctch[t-1];
	//	B[t] = B[t-1] + r*B[t-1]*(1-exp(log(B[t-1]) - log(K))) - ctch[t-1];
	//	logMean[t] = log(q) + log(B[t]);
	//}
	////
	//cpue ~ lognormal(logMean, s);
	////priors
	//q ~ beta(1, 5);
	//r ~ beta(3, 5);
	////logK ~ cauchy(0, 8);
	//K ~ cauchy(0, 2900);
	//s ~ cauchy(0, 1);
}
