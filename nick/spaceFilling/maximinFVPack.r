#
library(plgp)
library(boot)
library(mlegp)
library(pracma)

#
#FUNCTIONS
#

#
maximin = function(n, m, TT=100000) {
        #n : number of points to draw
        #m : dimension of space to draw from
        #
        X = matrix(runif(n*m), ncol=m)
        d = distance(X)
        d = d[upper.tri(d)]
        #
        md = min(d)
        for(t in 1:TT) {
                #propse a row
                row = sample(1:n, 1)
                xold = X[row,]
                #assume we accept
                X[row,] = runif(m)
                d = distance(X) #new distance
                d = d[upper.tri(d)]
                mdprime = min(d)
                if( mdprime > md ){ md=mdprime ## accept 
                } else { X[row,] <- xold } ## reject
        }
        return(X)
}

#
maximinAdd = function(n, m, TT=100000, XStart=NULL) {
        #n : number of points to draw
        #m : dimension of space to draw from
        #
        #value :  only retunrs new additions
        #
        X = matrix(runif(n*m), ncol=m)
        d = distance(X)
        d = d[upper.tri(d)]
        #
        md = min(d)
        #handle initiated design
        if(!is.null(XStart)){
                md2 = min(distance(X, XStart))
                if(md2 < md){ md=md2 }
        }
        for(t in 1:TT){
                #propse a row
                row = sample(1:n, 1)
                xold = X[row,]
                #assume we accept
                X[row,] = runif(m)
                d = distance(X) #new distance
                d = d[upper.tri(d)]
                mdprime = min(d)
                #handle initial design (make min dist smaller if neccesary with DS
                if(!is.null(XStart)){
                        mdprime2 = min(distance(X, XStart))
                        if(mdprime2 < mdprime) mdprime=mdprime2
                }
                #
                if( mdprime > md ){ md=mdprime ## accept 
                } else{ X[row,]=xold } ## reject
        }
        #return only the augmenting points
        return(X)
}


#
maximinFV = function(n, m, f, Xrange, Yrange, TT=100000, XStart=NULL, f_se=F) {
        #n : number of points to draw
        #m : dimension of space to draw from
        #f : function to spread over
        #Xrange : limits in X
        # I need to figure out how to actually place the limits on the space of y=f(x)'s
        #Yrange : limits in f(X)
        #TT     : samples to compute in maximin
        #XStart : a possible starting design to sequentially improve
        #f_se   : 
        #       if FALSE: maximin of distance of E(f)
        #       if TRUE: maximin of E(distance)
        #
        #value :  only retunrs new additions

        #
        Xrange = matrix(Xrange, nrow=m)
        X = matrix(rep(NA, n*m), ncol=m)
        #
        jobsNotDone = rep(T, n)
        while( any(jobsNotDone) ){
                #
                for(i in (1:n)[jobsNotDone]){
                        #
                        X[i,]=runif(m, Xrange[,1], Xrange[,2])
                }
                FV = rep(0, n)
                FF = f(X)
                if( f_se ){
                        FV = as.numeric(FF$se.fit^2)
                        FF = FF$fit
                }
                jobsNotDone = !(FF>Yrange[1] & FF<Yrange[2])
        }
        # 
        d = t(t(distance(FF)+FV)+FV)
        #
        md = min(d)
        #handle initiated design
        if(!is.null(XStart)){
                #
                FVStart = 0
                FFStart = f(XStart)
                if( f_se ){
                        #
                        FVStart = as.numeric(FFStart$se.fit^2)
                        FFStart = FFStart$fit
                }
                d = t(t(distance(FF, FFStart)+FV)+FVStart)
                md2 = min(d)
                if(md2 < md){ md=md2 }
        }
        for(t in 1:TT){
                #propose a row
                row = sample(1:n, 1)
                xold = X[row,]
                #
                jobNotDone = T
                while( jobNotDone ){
                        #assume we accept 
                        X[row,]=runif(m, Xrange[,1], Xrange[,2])
                        #
                        FV = rep(0, n)
                        FF = f(X)
                        if( f_se ){
                                FV = as.numeric(FF$se.fit^2)
                                FF = FF$fit
                        }
                        jobNotDone = !(FF[row]>Yrange[1] & FF[row]<Yrange[2])
                }
                #d = distance(FF) #new distance
                d = t(t(distance(FF)+FV)+FV)
                d = d[upper.tri(d)]
                mdprime = min(d)
                #handle initial design (make min dist smaller if neccesary with DStart included)
                if(!is.null(XStart)){
                        d = t(t(distance(FF, FFStart)+FV)+FVStart)
                        mdprime2 = min(d)
                        #mdprime2 = min(distance(FF, FFStart))
                        if(mdprime2 < mdprime) mdprime=mdprime2
                }
                #
                if( mdprime > md ){ md=mdprime ## accept 
                } else{ X[row,]=xold } ## reject
        }
        #return only the augmenting points
        return(X)
}

#
#2-D space filling( 
#gamma -> zeta
#(gamma, F*, M) -> alpha -> beta
#can pick F^*/M can pick gamma

#f: (gamma, alpha) -> (xi, zeta) <- fill this space
#given gamma alpha is known perfectly for a given xi
#g

#
getZeta = function(gamma, ff, M){
        (ff/(M+ff))/( gamma*((ff+M)/(M))^(1/gamma) - gamma + (ff/(M+ff)) )
}
getZeta = Vectorize(getZeta, "gamma")

#
f = function(gamma, ff, M){
	c(ff/M, getZeta(gamma, ff, M))
}

#
M = 0.2
xi = 2.5
ff = xi*M

#
isSE = T
#f = function(x) inv.logit(x) #x^3 #
gamLimits = c(-1, 1)
zetaLimits = getZeta(gamLimits, ff, M) #c(0, 1)
XStar = seq(gamLimits[1], gamLimits[2], length.out=5*(gamLimits[2]-gamLimits[1]))

#f will be made of GP model prediction
nStart = 10
#X = maximinF(nStart, 1, f, Xlimits, Ylimits, TT=100000, XStart=NULL)
gam = c(gamLimits[1], runif(nStart, gamLimits[1], gamLimits[2]), gamLimits[2])
zeta = getZeta(gam, ff, M)
#
Xall = t(t(gam))
Yall = t(t(zeta))

#
M = 50-nStart
for(i in 1:M){
        #FIT and update stat model
        gpFit = mlegp(Xall, Yall, constantMean=F, verbose=0)
        fHat = function(x, d=1){
                #
                x = matrix(x, ncol=d)
                pred = predict(gpFit, x, se.fit=isSE)
                #
                return( pred )
        }
#
#        #plot new predicitons
#        FF = fHat(XStar)
#        if(isSE){ FF=FF$fit }
#        lines(XStar, FF, col=i)
#
#        #find new points to evaluate
#        X = maximinFV(1, 1, fHat, Xlimits, Ylimits, TT=25000, XStart=Xall, f_se=isSE)
#        y = f(X)
#        #add nnew evaluations to all evaluations
#        Xall = rbind(Xall, X)
#        Yall = rbind(Yall, y)
#
#        ##plot new evaluations   
#        #text(X, y, label=as.character(nStart+i), col=i)
#        #text(X, rep(Ylimits[1], length(y)), label=as.character(nStart+i), col=i)
#        #text(rep(Xlimits[1], length(X)), y, label=as.character(nStart+i), col=i)
#}

