data{
        int<lower=0> N;
        real<lower=0> cpue[N];
        int<lower=0> C[N];
}

parameters{
        real<lower=0> q;
        real<lower=0> alpha;
        real<lower=0> beta;
	real<lower=0, upper=1> sigma;
        real<lower=0> s;
}

model{
        //
	int k;
	int h;
	//
	real wk;
	real sh;
        real gamma;
        real kappa;
	//
	real W0;
	real rh;
	real R0; 
        real B0;
        real C0;
	//
        real B[N];
        real W[N];
        real R[N];
        real logI[N];
	
	//
	k = 2;
	h = 0;
	//
	wk = 0.001;
	sh = sigma*(1-h);
	gamma = -1;
	kappa = 0;	
	// 
        W0 = ((1-sh)*wk + (1-kappa)*sh*wk)/(1-kappa*sigma);
        rh = (wk/W0)*(1-sh)/(1-h);
        R0 = (1-(rh/alpha)^gamma)*rh/(beta*gamma); 	//#(1-exp(gamma*(log(rh)-log(alpha))))*(rh/(beta*gamma)) #
        B0 = (W0/wk)*(R0/(1-sh)); 			//#(W0*Rh)/(wk*#b#10000;
        C0 = h*B0;
	
	//
	R[1] = R0; 
        B[1] = R[1] + (((1-kappa)*wk + kappa*W0)/W0) * sigma*(B0-C0); //#exp(log( (1-kappa)*wI+kappa*W0 ) + lo
        W[1] = B[1]/( R[1]/wk + sigma*(B0-C0)/W0 );
        logI[1] = log(q) + log(B[1] - C[1]);
        for(t in 2:N){
                if(t>k){ 
			R[t] = alpha*(B[t-k]-C[t-k])*(1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma);
                }else{   
			R[t] = R0; 
		}
                B[t] = R[t] + (((1-kappa)*wk + kappa*W[t-1])/W[t-1]) * sigma*(B[t-1]-C[t-1]); //#exp(log( (1-k
                W[t] = B[t]/( R[t]/wk + sigma*(B[t-1]-C[t-1])/W[t-1] );
                logI[t] = log(q) + log(B[t] - C[t]);
        }
	
	//
	cpue ~ lognormal(logI, s);
	//
        q ~ cauchy(0, 1e-3);
        alpha ~ cauchy(0, 1e1);
        beta ~ cauchy(0, 1e-2);
	sigma ~ beta(5, 1);
        s ~ cauchy(0, 0.5);

}
