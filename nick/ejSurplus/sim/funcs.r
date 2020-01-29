#
abGet = function(hs, Cs, delta, gamma, wk, kappa){
	#hs:
	#Cs:
	#delta:
	#gamma:
	#wk:
	#kappa:
	#
	sh = (1-delta)*(1-hs)
	Ws = ((1-sh)*wk + (1-kappa)*sh*wk) / (1-(kappa*sh))
	Qs = ((1-(wk/Ws))*(1-sh)) / (1-(kappa*sh))
	#
	alpha = (wk/Ws) * ((1-sh)/(1-hs)) * (1+((gamma*(1+Qs)*hs)/(1-sh)))^(1/gamma)
	beta  = ((1+Qs)*hs^2) / (Cs*(1-hs)*(1-sh+gamma*hs*(1+Qs)))
	#
	return( data.frame(a=alpha, b=beta) )
}

#
sGet = function(h, delta){
	sh = (1-delta)*(1-h)
	return( sh )	
}

#
WGet = function(s, wk, kappa){
	W = ((1-s)*wk + (1-kappa)*s*wk) / (1-(kappa*s))
	return( W )
}

#
rGet = function(h, s, W, wk){
	r = (wk/W)*((1-s)/(1-h))
	return( r )
}

#
RGet = function(alpha, beta, gamma, r){
	R = (1-(r/alpha)^gamma) * r/(beta*gamma)
	return( R )
}

#
BGet = function(s, W, R, wk){
	B = (W/wk)*(R/(1-s))
	return( B )
}

#
f = function(x, pars){
        #
        Cs    = x[1]
        gamma = x[2]
        #
        kappa = pars$kappa
        wk    = pars$wk
        wI    = pars$wI
        Bs    = pars$Bs
        B0    = pars$B0
        M     = pars$M
        hs    = pars$hs
        #
        delta = 1-exp(-M)
        #
        ab    = abGet(hs, Cs, delta, gamma, wk, kappa)
        alpha = ab$a
        beta  = ab$b
        #
        h  = 0
        s0 = sGet(h, delta)
        W0 = WGet(s0, wk, kappa)
        r0 = rGet(h, s0, W0, wk)
        R0 = RGet(alpha, beta, gamma, r0)
        B0E = BGet(s0, W0, R0, wk)-B0
        #
        ss = sGet(hs, delta)
        Ws = WGet(ss, wk, kappa)
        rs = rGet(hs, ss, Ws, wk)
        Rs = RGet(alpha, beta, gamma, rs)
        BsE = BGet(ss, Ws, Rs, wk)-Bs
        #
        return( c(BsE, B0E) )
}

#
model = function(q, alpha, beta, delta, kappa, gamma, Y, C){
        k = 2
        wI = 0.001
        wk = (1-kappa)*wI + kappa*wI#w_(a-1)
        #
        h  = 0
        s0 = sGet(h, delta)
        W0 = WGet(s0, wk, kappa)
        r0 = rGet(h, s0, W0, wk)
        R0 = RGet(alpha, beta, gamma, r0)
        B0 = BGet(s0, W0, R0, wk)
        C0 = h*B0
        #
        B = matrix(NA, Y, 1)
        W = matrix(NA, Y, 1)
        R = matrix(NA, Y, 1)
        logI = matrix(NA, Y, 1)
        ##
	#print(alpha)
	#print(beta)
	#print(kappa)
	#print(delta)
	#print(wI)
	#print(wk)
	##
        flag = (alpha*( 1 + ( (1-kappa)*(1-delta) )/( delta+(1-kappa)*(1-delta) ) * ((wI/wk) - 1) ))>(delta)
        if( !flag ){ return(NA) }
        #
        R[1] = R0 #alpha*(B0-C0)*(1-beta*gamma*(B0-C0))^(1/gamma)		
        B[1] = R[1] + (((1-kappa)*wI + kappa*W0)/W0) * (1-delta)*(B0-C0) #exp(log( (1-kappa)*wI+kappa*W0 ) + log(sigma) + log(B0) - log(W0	
        W[1] = B[1]/( R[1]/wk + (1-delta)*(B0-C0)/W0 );	 
	logI[1] = log(q) + log(B[1] - C[1]);
	#	
        for(t in 2:Y){	
                if(t>k){ R[t] = alpha*(B[t-k]-C[t-k])*(1-beta*gamma*(B[t-k]-C[t-k]))^(1/gamma)
                }else{   R[t] = R0 }
                B[t] = R[t] + (((1-kappa)*wI + kappa*W[t-1])/W[t-1]) * (1-delta)*(B[t-1]-C[t-1]) #exp(log( (1-kappa)*
                W[t] = B[t]/( R[t]/wk + (1-delta)*(B[t-1]-C[t-1])/W[t-1] );
                logI[t] = log(q) + log(B[t] - C[t]);
        }
        #
        out = list(
                R=R,
                B=B,
                W=W,
                logI=logI
        )
        #
        return( out )
}



