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
				"iterate", "optimize", "model",	"prior", "like",
				"plotQuan", "plotMean", "plotBand", "plotRS", 
				"quan", "N0Funk", "B0Funk", "save", "load", "printer"
			))
		},
		plotQuan = plotQuan,
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


##
##TEST
##
#
##
#w = function(a, wi, k){ wi*(1-exp(-k*a)) }
#
##
#f = function(t, Y, lalpha, lbeta, gamma, a0, wi, k, catch, B0){
#        #linearly interpolate catches
#        ft = floor(t)
#        q  = (t-ft)
#        Cl = catch[ft]
#        Cu = catch[ft+1]
#        C = q*Cu + (1-q)*Cl
#        if(q==0){ C=Cl }
#	#
#        N = Y[1]
#        B = Y[2]
#        #
#        if( (t-a0)<1){
#                Blag = B0
#        }else{
#                Blag = lagvalue(t-a0)[2]
#        }
#        #
#	#R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
#        alpha = exp(lalpha)
#	beta = exp(lbeta)
#	R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
#	#print(C)
#        out = c(N=NA, B=NA)
#        out[1] = R - (M+C)*N
#        out[2] = wi*(1-exp(-k*a0))*R + k*(wi*N-B) - (M+C)*B
#        #
#        return( list(out) )
#}
#
##
#g = function(t, Y, p){ #alpha, beta, gamma, a0, wi, k, catch){
#        #linearly interpolate catches
#        ft = floor(t)
#        q  = (t-ft)
#        Cl = catch[ft]
#        Cu = catch[ft+1]
#        C = q*Cu + (1-q)*Cl
#        if(q==0){ C=Cl }
#	#
#        N = Y[1]
#        B = Y[2]
#        #
#        if( (t-a0)<1){
#                Blag = P0
#        }else{
#                Blag = lagvalue(t-a0)[2]
#        }
#        #
#	#R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
#        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
#	#
#        out = c(N=NA, B=NA)
#        out[1] = R - (M+C)*N
#        out[2] = wi*(1-exp(-k*a0))*R + k*(wi*N-B) - (M+C)*B
#        #
#        return( list(out) )
#}
#
##
#Fs = c(seq(0.2, 2, length.out=15), rev(seq(0.1, 2, length.out=15)), rep(0.1, 15))
#catch = Fs
#
##
#M = 0.2
#k = 0.2
#
##
#wi = 1
#a0 = 2
#TT = 45
#
##
#alpha = 5
#beta = 1
#gamma = -0.99
##gamma=0 is a problem
##gamma=-1 is BH in limit; a problem otherwise
#
##
#Rtil = alpha/(beta*(1+gamma)) * (1-gamma/(1+gamma))^(1/gamma)
#N0 = Rtil/(M)#+FF)
#P0 = (w(a0, wi, k)*Rtil + k*wi*N0)/(k+M)#+FF)
#
##
#dOut = dede(c(N0, P0), 1:TT, g, NULL, method="lsode")
#
##
#dat = ddModel$new( derivs=f,
#        N0Funk=function(lalpha, lbeta, gamma){
#                #
#                alpha = exp(lalpha)
#                beta  = exp(lbeta)
#                #
#                ( alpha/(beta*(1+gamma)) * (1-gamma/(1+gamma))^(1/gamma) )/M
#        },
#        B0Funk=function(lalpha, lbeta, gamma, wi, k){
#                #
#                alpha = exp(lalpha)
#                beta = exp(lbeta)
#                #
#                ( wi*(1-exp(-k*a0))*(alpha/(beta*(1+gamma))*(1-gamma/(1+gamma))^(1/gamma)) +
#                  k*wi*(alpha/(beta*(1+gamma))*(1-gamma/(1+gamma))^(1/gamma))/M
#                )/(k+M)
#        },
#        time=1:TT, catch=Fs, a0=a0, M=M, wi=wi, k=k,    #constants
#        lalpha=log(alpha), lbeta=log(beta), gamma=0.2, #parameters
#        lq=log(0.00049), lsdo=log(0.01160256)           #nuisance parameters
#)
#dat$iterate()
#
##
#test = ddModel$new( derivs=f, 
#	N0Funk=function(lalpha, lbeta, gamma){
#		#
#		alpha = exp(lalpha)
#		beta  = exp(lbeta) 
#		#
#		( alpha/(beta*(1+gamma)) * (1-gamma/(1+gamma))^(1/gamma) )/M 
#	},
#        B0Funk=function(lalpha, lbeta, gamma, wi, k){ 
#		#
#		alpha = exp(lalpha)
#		beta = exp(lbeta)
#		#
#		( wi*(1-exp(-k*a0))*(alpha/(beta*(1+gamma))*(1-gamma/(1+gamma))^(1/gamma)) + 
#		  k*wi*(alpha/(beta*(1+gamma))*(1-gamma/(1+gamma))^(1/gamma))/M
#		)/(k+M)
#	},
#	time=1:TT, catch=Fs, a0=a0, M=M, wi=wi, k=k, 	#constants
#        lalpha=log(alpha), lbeta=log(beta), gamma=gamma,#parameters
#        lq=log(0.00049), lsdo=log(0.01160256)  		#nuisance parameters
#)
#test$iterate()
#
###
##cpue  = c(1.78, 1.31, 0.91, 0.96, 0.88, 0.90, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49, 0.43, 0.40, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63)
##catch = c(94, 212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254, 170, 97, 91, 177, 216, 229, 211, 231, 223)
##
##test$time = 1:length(cpue)
##test$catch = catch
##test$iterate()
#
##
#library(pracma)
#
##
#d = exp( rnorm(TT, dat$lq+log(dat$B), exp(dat$lsdo)) ) #exp(log(dat$B)+dat$lq+rlnorm(TT)
#
##
#test$optimize( d,
#	c('lsdo', 'lalpha', 'lbeta'),	
#	lower = c(-10, log(eps()), log(eps())),
#	upper = c(10, log(10), log(10)),
#	gaBoost = list(run=10, parallel=T, popSize=5*10^3),
#	fitQ = F,
#	cov = T
#)
#test$printSelf()
#
##
#plot(d)
#test$plotMean(add=T)
#test$plotBand()
##
#dat$plotMean(add=T, col='red')
#
###
##test$optimize( d,
##	c('lsdo', 'lalpha', 'lbeta', 'gamma'),	
##	lower = c(-10, log(eps()), log(eps()), -1),
##	upper = c(10, log(10), log(10), 2),
##	gaBoost = list(run=10, parallel=T, popSize=5*10^3),
##	fitQ = F,
##	cov = T
##)
##
####
###test$optimize( d,
###        c('lsdo', 'lalpha', 'lbeta'),
###        lower = c(-10, log(eps()), log(eps())),
###        upper = c(10, log(10), log(10)),
###        #gaBoost = list(run=10, parallel=T, popSize=5*10^3),
###        fitQ = T,
###        cov = T
###)
##
###
##test$printSelf()
##test$plotMean(add=T, col='blue')
#
##
#dev.new()
#test$plotQuan(function(N){N})
#dev.new()
#test$plotQuan(function(N, B){B/N})
