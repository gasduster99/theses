data{
	int<lower=0> N;
	real<lower=0> cpue[N];
	//'catch' is an internal variable to stan; no variable names are allowed to be call 'catch'
	real<lower=0> ctch[N];
}

parameters{
	//real<lower=0, upper=1> q;
        //real<lower=0, upper=1> r;
        //real<lower=0> K;
        real<lower=0> s;
        real<lower=0> B0;
}

model{
	//q prior
	//r prior
	//K prior
	//s prior
	//B0 prior

	//log(cpue) ~ normal(log(I), s);
	cpue ~ lognormal(log(B0), s);
}
