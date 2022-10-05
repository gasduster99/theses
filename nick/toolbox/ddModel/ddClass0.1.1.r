#
suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(deSolve, quietly=TRUE) ))
#
source('modelFunk.r')
source('outputFunk.r')
source('classifyFunk.r')
source('parSelfParFunk.r')
source('optimizeFunk.r')

#
#CLASS
#

#
ddModel = R6Class("DDModel", lock_objects=FALSE,
	#
	public = list(
		#pop
		N  = NA,
		N0 = NA,
		B  = NA,
		B0 = NA,	
		N0Funk = NA,
		B0Funk = NA,
		#time
		time = NA,
		#model
                lsdo = NA, 
		lq = 0, #log proportionality constant between cpue and N
		model = list(observation="LN"),
                prior = list(),
		#functions
		#computational
		ODE_method = 'lsode',
		OPT_method = 'L-BFGS-B',
		
		#
		initialize = function( 	N0   = NA,
					B0   = NA,
					time = NA,
					derivs = NA,
					N0Funk = NA,
					B0Funk = NA,
					...
		){	
			#misc variable digestion
                        misc = list(...)
                        miscNames = names(misc)
                        for(m in miscNames){
                                eval(parse( text=sprintf("self$%s=misc[[m]]", m) ))
                        }
			
			#dNdt
			stopifnot(is.function(derivs))
			private$dNdt_classify(derivs)	
			
			#N0
			if( is.function(N0Funk) ){ 	
				private$N0_classify(N0Funk)
				N0 = self$N0Funk()
			}
			#B0
			if( is.function(B0Funk) ){ 	
				private$N0_classify(B0Funk)
				B0 = self$B0Funk()
			}		
	
			#preallocate N
			self$N0 = N0
			self$B0 = B0
			self$time = time
			self$N = matrix(NA, nrow=length(time), ncol=1)
			self$B = matrix(NA, nrow=length(time), ncol=1)
			rownames(self$N) = sprintf("TIME %d", time)
			rownames(self$B) = sprintf("TIME %d", time)
		},
		
		#
		iterate = function(method=self$ODE_method){
			#method	: optionally update the ODE method to be handed to the ode function (default 'rk4')
			#
			#value	: nothing returned, but self$N is updated with current values

			#prechecking 
			
			#digest and possibly change ode method	
			self$ODE_method = method	
	
			#N0 & B0
                        if( is.function(self$N0Funk) ){ self$N0=self$N0Funk() }
			if( is.function(self$B0Funk) ){ self$B0=self$B0Funk() }
	
			#last minute allocation and variable updates
			#self$N0 = N0
			#self$B0 = B0
                        self$N = matrix(NA, nrow=length(self$time), ncol=1)
			self$B = matrix(NA, nrow=length(self$time), ncol=1)
                        rownames(self$N) = sprintf("TIME %d", self$time)
			rownames(self$B) = sprintf("TIME %d", self$time)		

			#solve 
        		#capture.output( self$N <- ode(self$N0, self$time, private$dNdt, parms=NULL, method=method)[,2], file="/dev/null" )
			capture.output( out <- dede(c(self$N0, self$B0), self$time, private$dNdt, parms=NULL, method=method), file="/dev/null" )
			self$N = out[,2]
			self$B = out[,3]
		},
		
		#
                optimize = optimize,
		
		#
		printer   = printSelf, 
		printSelf = function(ins=c()){
			self$printer(ins, outs=c(
				"iterate", "optimize", "model",	"prior", 
				"plotMean", "plotBand", "plotRS", "N0Funk", 
				"save", "load", "printer"
			))
		},
		plotMean = plotMean,
		plotBand = plotBand,
		plotRS 	 = plotRS,
		#
		save = function(fileName){ saveRDS(self, file=fileName) },
		load = function(fileName){ readRDS(fileName) },
		#
		like = function(data){ sum(private$dLikes[[self$model$observation]](self, data)) }
	),
	
	#
	private = list(
		##
		#dNdt = NA,
		
		#
		selfToPar = selfToPar,	
		parToSelf = parToSelf, #NOTE: parValues should be passed with names
		
		#
		dNdt_classify = dNdt_classify,
		N0_classify = N0_classify,
		classify = classify,
		
		#
		dLikes = dLikes,
		qLikes = qLikes
	)
)


#
#TEST
#

#
w = function(a, wi, k){ wi*(1-exp(-k*a)) }

#
f = function(t, Y, alpha, beta, gamma, a0, wi, k, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
	#
        N = Y[1]
        B = Y[2]
        #
        if( (t-a0)<1){
                Blag = 0
        }else{
                Blag = lagvalue(t-a0)[2]
        }
        #
	#R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
	#
	print(C)
        out = c(N=NA, B=NA)
        out[1] = R - (M+C)*N
        out[2] = wi*(1-exp(-k*a0))*R + k*(wi*N-B) - (M+C)*B
        #
        return( list(out) )
}

#
g = function(t, Y, p){ #alpha, beta, gamma, a0, wi, k, catch){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
	#
        N = Y[1]
        B = Y[2]
        #
        if( (t-a0)<1){
                Blag = 0
        }else{
                Blag = lagvalue(t-a0)[2]
        }
        #
	#R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
	#
        out = c(N=NA, B=NA)
        out[1] = R - (M+C)*N
        out[2] = wi*(1-exp(-k*a0))*R + k*(wi*N-B) - (M+C)*B
        #
        return( list(out) )
}

#
P0 = 10000
N0 = 100 #P0/1000
wi = P0/N0*10

#
a0 = 10
TT = 45
#
Fs = c(seq(0.2, 2, length.out=15), rev(seq(0.1, 2, length.out=15)), rep(0.1, 15))
M = 0.2
k = 0.2

#
catch = Fs
alpha=1; beta=1; gamma=-1
dOut = dede(c(N0, P0), 1:TT, g, NULL, method="lsode")

#
test = ddModel$new(
        derivs=f, N0=N0, B0=P0, #N0Funk=function(){100000}, P0Funk=function(){10000}, #model
        time=1:TT, catch=Fs, a0=a0, M=M, wi=wi, k=k, 	#constants
        alpha=alpha, beta=beta, gamma=gamma,    	#parameters
        lq=log(0.00049), lsdo=log(0.01160256)  		#nuisance parameters
)
