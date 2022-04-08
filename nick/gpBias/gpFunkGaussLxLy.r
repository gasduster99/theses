#
library(boot)
library(MASS)
library(pracma)
library(mvtnorm)
library(parallel)
library(extraDistr)

#
#Functions for fitting a GP
#

#
yStarOnGrid = function(axes, axeGrid, yGrid, p1=1, p2=1){
        #axes   : nx2 matrix with columns 
        #axeGrid: px2 matrix of predicted locations
        #yGrid  : p-vector of predictions over xGrid
        #

        #
        ax1 = round(axes[,'x^1y^0'], p1)
        ax2 = round(axes[,'x^0y^1'], p2)
        #
        yOut = matrix(NA, nrow=nrow(axes), ncol=1)
        for(i in 1:nrow(axes)){
                where = which(axeGrid[,1]==ax1[i] & axeGrid[,2]==ax2[i])
                yOut[i] = yGrid[where]
        }
        #
        return( yOut )
}

#
strongInv = function(X){
        out = tryCatch({
                #
                chol2inv(chol(X))
        #
        }, warning=function(war){
                #
                print(war)
                #
                return(war)
        #
        }, error=function(err){
                #
                print(err)
                writeLines('Trying SVD Moore-Penrose ginv')
                #       
                #SVD = svd(SIG)
                #ISIG = SVD$v%*%diag(1/SVD$d)%*%t(SVD$u) 
                return(ginv(X))
        })
        #
        return(out)
}

#
logDet = function(x){ determinant(x, log=T)$mod[1] }

##
#maternCor = function(t, v){
#        #d = t/phi
#	d = t
#        d[t == 0] = eps()
#        con = (2^(v-1))*gamma(v)
#        con = 1/con
#        out = (con*(d^v)*besselK(d, v))*(0<t & 0<v) + as.numeric(t==0)
#        #
#        return( out )
#}

#
powExpCor = function(t, phi, v=2){ exp(-abs(t/phi)^v)*(0<t & 0<v & v<=2) + as.numeric(t==0) }

#
logML = function(x, fb){
	#
	l1=x['l1']; l2=x['l2']; s2=x['s2']
	#     from like             from beta integral
	out = -0.5*logDet(fb$V) + -0.5*logDet(t(fb$LIX)%*%fb$LIX) + -0.5*fb$S2 + 
		s2Prior(s2) + l1Prior(l1) + l2Prior(l2) #+ thPrior(th) + nuPrior(nu)
	#
	return(out)
}

#
funnyBusiness = function(l1, l2, s2){	
	#
	SIG = s2*powExpCor(dMat1, l1)*powExpCor(dMat2, l2)	
	#
	L  = t(chol(SIG))
        LIy = forwardsolve(L, y)
        LIX = forwardsolve(L, X)
	#
	bHat = .lm.fit(LIX, LIy)$coef
	out = list(
                V  = SIG,
                LIX  = LIX,
                bHat = bHat,
                S2   = t(LIy-LIX%*%bHat)%*%(LIy-LIX%*%bHat)
        )
        return( out )
}

#
rectDist = function(x, y){
        mcmapply(function(y1){
                apply(x, 1, function(xi){
                        norm(y1-xi, '2')
                })
        }, y, mc.cores=detectCores())
}

#
l1a=2; l1b=1 #0.1
l1Prior = function(l1, log=T){
        #
        out = dinvgamma(l1, l1a, l1b, log=T)
        if( !log ){
                out = dinvgamma(l1, l1a, l1b, log=F)
        }

        #
        return( out )
}

#
l2a=2; l2b=1
l2Prior = function(l2, log=T){
        #
        out = dinvgamma(l2, l2a, l2b, log=T)
        if( !log ){
                out = dinvgamma(l2, l2a, l2b, log=F)
        }

        #
        return( out )
}

#
s2a=2; s2b=0.1
s2Prior = function(s2, log=T){
	#
        out = dinvgamma(s2, s2a, s2b, log=T)
        if( !log ){
                out = dinvgamma(s2, s2a, s2b, log=F)
        }

        #
        return( out )
}

#
registerData = function(y, X, axes, Tg, dMat1=dist(axes[,1]), dMat2=dist(axes[,2]), l1a=.GlobalEnv$l1a, l1b=.GlobalEnv$l1b, l2a=.GlobalEnv$l2a, l2b=.GlobalEnv$l2b, s2a=.GlobalEnv$s2a, s2b=.GlobalEnv$s2b){ #nua=.GlobalEnv$nua, nub=.GlobalEnv$nub, tha=.GlobalEnv$tha, thb=.GlobalEnv$thb
	#Hack the system and force the global envirnment to contain these variables for easy passing to functions above.
	#dangerous if you don't know its happening, but your reading this and I wrote it so we all know its happening, 
	#sooo... WE'RE GOOD. Get over it!
	.GlobalEnv$y = y
	.GlobalEnv$X = X
	.GlobalEnv$axes = axes
	.GlobalEnv$Tg = Tg
	#
	.GlobalEnv$dMat1 = as.matrix(dMat1)
	.GlobalEnv$dMat2 = as.matrix(dMat2)
	#
	.GlobalEnv$l1a=l1a; .GlobalEnv$l1b=l1b
	.GlobalEnv$l2a=l2a; .GlobalEnv$l2b=l2b	
	.GlobalEnv$s2a=s2a; .GlobalEnv$s2b=s2b
}

#
gpMAP = function( par, lower=rep(eps(), 3), upper=c(Inf, Inf, Inf), hessian=F, psiSample=F ){	
	#
	realPar = par
	realPar['l1'] = log(realPar['l1'])
	realPar['l2'] = log(realPar['l2'])	
	realPar['s2'] = log(realPar['s2'])
	#
	opt = optim(realPar,
        	function(lp){
			p = c(exp(lp['l1']), exp(lp['l2']), exp(lp['s2']))
			fb = funnyBusiness(p['l1'], p['l2'], p['s2']) 
			return(-logML(p, fb))
        	},
        	control=list('maxit'=10^6),
        	hessian=(hessian | psiSample!=F), 
        	method="L-BFGS-B",
        	lower=log(lower),
        	upper=log(upper)
	)
	#
	outPar = opt$par
	outPar['l1'] = exp(outPar['l1'])
	outPar['l2'] = exp(outPar['l2']) 	
	outPar['s2'] = exp(outPar['s2'])
	#
	fb = funnyBusiness(outPar['l1'], outPar['l2'], outPar['s2'])
	out = list(beta=fb$bHat, psi=outPar)
	#
	if( hessian | psiSample!=F ){ 
		#
		psiCov = chol2inv(chol(opt$hessian))
		betaCov = chol2inv(chol(t(fb$LIX)%*%fb$LIX))
		#
		howMany = 10^4
		if(typeof(psiSample)=="double"){ howMany=psiSample }
		#
		sam = rmvnorm(howMany, opt$par, psiCov)
		sam[,'l1'] = exp(sam[,'l1'])
        	sam[,'l2'] = exp(sam[,'l2'])
        	sam[,'s2'] = exp(sam[,'s2'])		
		#
		out$betaCov = betaCov
		out$psiCov = cov(sam)
		#
		if( psiSample!=F ){ out$psiSam=sam }
	}
	
	#
	return(out)
}

#
gpPredict = function(XStar, axeStar, gpFit, asMat=T){
	#
	XStar = as.matrix(XStar)
	SIGStar = powExpCor(rectDist(as.matrix(axeStar[,1]), as.matrix(axes[,1])), gpFit$psi['l1'])
	SIGStar = SIGStar * powExpCor(rectDist(as.matrix(axeStar[,2]), as.matrix(axes[,2])), gpFit$psi['l2'])
	#
	SIG = powExpCor(dMat1, gpFit$psi['l1'])*powExpCor(dMat2, gpFit$psi['l2'])
	#
	out = XStar%*%gpFit$beta + SIGStar%*%strongInv(SIG)%*%(y-X%*%gpFit$beta)
	#
	if( asMat ){
		out = matrix(out, nrow=length(unique(axeStar[,1])), ncol=length(unique(axeStar[,2]))) 
	}
	#
	return( out )
}

##
##Test
##
#
###
##y = 1:50
##X = matrix(c(rep(1,50), 1:50, seq(0.1, 100, 2.1)), ncol=3)
##axes = X[,c(2,3)]
###
##Tg = diag(1:50/10)
#
##
#library(geoR)
#library(stringr)
##
#dat = read.csv('~/Documents/school/ucscGrad/stat226/hw4/D.csv')[,-1]
#set.seed(2)
#dat[,1:2] = jitterDupCoords(cbind(dat$long, dat$lat), max=0.1) #, max=0.00000000001)
##
#deg = 2
#pol = poly(dat$long, dat$lat, degree=deg, raw=T)
#X = cbind(1, pol)
#polNames = str_replace( sprintf('x^%s', colnames(pol)), '\\.', 'y^' )
#colnames(X) = c('x^0y^0', polNames)
##
#y = logit(dat$moist)
#n = length(y)
#k = ncol(X)
#axes = X[,c('x^1y^0', 'x^0y^1')]
#
##
#Tg = diag(n)*0.0123/2
#
##
#par=c(l1=0.5, l2=0.5, s2=0.28)
#registerData(y, X, axes, Tg)
#
##
#gpFit = gpMAP(par, hessian=T, psiSample=F)
#
##
#lonGrid = seq(-105, -89, 0.1)
#latGrid = seq(42, 50, 0.1)
#XStar = expand.grid(lonGrid, latGrid)
##
#deg = 2
#pol = poly(XStar[,1], XStar[,2], degree=deg, raw=T)
##
#XStar = cbind(1, pol)
#polNames = str_replace( sprintf('x^%s', colnames(pol)), '\\.', 'y^' )
#colnames(XStar) = c('x^0y^0', polNames)
#axeStar = XStar[,c('x^1y^0', 'x^0y^1')]
#
##
#gpPred = gpPredict(XStar, axeStar, gpFit)
#gpMat = matrix(gpPred, nrow=length(lonGrid), ncol=length(latGrid))
#
##
#image(lonGrid, latGrid, gpMat)






#
#JUNK
#


#fb = funnyBusiness(par['l1'], par['l2'], par['th'], par['nu'], par['s2'])
#ml = logML(par, fb)
#opt = optim(par, 
#	function(par){
#		fb = funnyBusiness(par['l1'], par['l2'], par['th'], par['nu'], par['s2'])
#		return(-logML(par, fb))
#	}, 
#	control=list('maxit'=10^6), 
#	#hessian=T, 
#	method="L-BFGS-B", 
#	lower=c(rep(eps(), 2),0,rep(eps(), 2)), 
#	upper=c(1, 1, pi/2, 3, 1)
#)
#fb = funnyBusiness(opt$par['l1'], opt$par['l2'], opt$par['th'], opt$par['nu'], opt$par['s2'])


##
#funnyBusiness = function(gamma, phi, v){
#        #
#        SIG  = maternCor(dis, phi, v)/gamma + diag(n)
#        #
#        L  = t(chol(SIG))
#        LIy = forwardsolve(L, y)
#        LIX = forwardsolve(L, X)
#        #
#        bHat = .lm.fit(LIX, LIy)$coef
#        out = list(
#                SIG  = SIG,
#                LIX  = LIX,
#                bHat = bHat,
#                S2   = t(LIy-LIX%*%bHat)%*%(LIy-LIX%*%bHat)
#        )
#        return( out )
#}

##
#logML = function(x, fb){
#        #
#        gamma = x['gamma']
#        phi   = x['phi']
#        v     = x['v']
#
#        #      from like             from beta integral                from t2 integral jac                   jacobian
#        out  = -0.5*logDet(fb$SIG) + -0.5*logDet(t(fb$LIX)%*%fb$LIX) + -((n-k)/2 + at + -1)*log(fb$S2/2+bt) + -gamma
#
#        #
#        return( out )
#}

##
#nua=2; nub=3
#nuPrior = function(nu, log=T){
#        #
#        out = dinvgamma(nu, nua, nub, log=T)
#        if( !log ){
#                out = dinvgamma(nu, nua, nub, log=F)
#        }
#
#        #
#        return( out )
#}
#
##
#tha=1; thb=3
#thPrior = function(th, log=T){
#	#
#        out = dbeta(th/(pi/2), tha, thb, log=T) - log(pi/2)
#        if( !log ){
#                out = dbeta(th/(pi/2), tha, thb, log=F)/(pi/2)
#        }
#
#        #
#        return( out )
#}
#
##NOTE:rewrite/consider
#anisoNorm = function(AX1, AX2, l1, l2, th){
#	#AX: an mx2 matrix (this function is designed to be more efficient with nrow(AX1)>nrow(AX2))
#	#AX: an nx2 matrix
#	#l1: a scalar \in R+ length scale in the 1st dimention
#	#l2: a scalar \in R+ length scale in the 2nd dimention
#	#th: a scalar \in (0, pi/2) rotation
#	
#	#
#	P = matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)
#	LI = diag(1/c(l1, l2))
#	KI = P%*%LI%*%t(P)
#	##serial
#	#mapply(function(y1, y2){
#	#	apply(AX, 1, function(xi){
#        #       		sqrt(t(c(y1, y2)-xi)%*%KI%*%(c(y1, y2)-xi))
#        #        })
#	#}, AX[,1], AX[,2])
#	#parallel
#	mcmapply(function(y1, y2){
#		apply(AX1, 1, function(xi){
#               		sqrt(t(c(y1, y2)-xi)%*%KI%*%(c(y1, y2)-xi))
#                })
#	}, AX2[,1], AX2[,2], mc.cores=detectCores()-1) #mc.preschedule=F,
#}

