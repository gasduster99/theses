#
suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(mvtnorm, quietly=TRUE) ))

#
source("outputFunk.r")
source("parSelfParFunk.r")

#
#CLASS
#

#
gpModel = R6Class("GPModel", lock_objects=FALSE,
        #
        public = list(
		#
		X = NA,
		Y = NA,
		g = 0,
		#
		initialize = function(
			S2 = NA,
			...
		){
			#misc variable digestion
                        misc = list(...)
                        miscNames = names(misc)
                        for(m in miscNames){
                                eval(parse( text=sprintf("self$%s=misc[[m]]", m) ))
                        }
			
			#
			self$yRes = self$Y-predict(self$lm)
			
			#
			if(!'obsV'%in%miscNames){ self$obsV=rep(0, length(self$Y)) }
			
			#Covariance function
                        stopifnot(is.function(S2))
                        private$classify(S2, 2)
		},
			
		#
		fit = function(	parNames,
                        	lower, upper,
                        	cov=F,
                        	control=list()
		){
			#
			loglikeL = function(par){
			        #
				private$parToSelf(par)
				#
			        SIG2 = self$S2(self$X, self$X) + diag(self$obsV) + self$g*diag(length(self$Y))
				-dmvnorm(self$yRes, sigma=SIG2, log=T)
			}
			#
			out = optim(
                        	private$selfToPar(parNames),
                        	loglikeL,
                        	lower = lower,
                        	upper = upper,
                        	hessian = cov,
                        	method  = 'L-BFGS-B',
                        	control = control
                	)
			#
			self$aic = 2*(out$value + length(parNames))
			#
                	if( cov ){
                	        #
                	        self$rsCov = chol2inv(chol(out$hessian))
                	        colnames(self$rsCov) = parNames
                	        rownames(self$rsCov) = parNames
                	}
			return( out )
		},
		
		#
		printer = printSelf,
                printSelf = function(ins=c()){
                        self$printer(ins, outs=c(
                                "S2", "fit", "save", "load", "lm", "printer"
                        ))
                },
		
		#
                save = function(fileName){ saveRDS(self, file=fileName) },
                load = function(fileName){ readRDS(fileName) }
	),
	#
	private = list(
		#
		classify = function(fun, numPars=1){
		        #
		
		        #
		        name = as.character(substitute(fun))
		        arg = formals(fun)
		        bdy = body(fun)[-1]
		        #get the names of arguments to replaces with self variables
			mainArgs = names(arg[1:(numPars)])
		        xtraArgs = names(arg[(numPars+1):length(arg)])
			#get arguments with defaults
			defArgs = arg[xtraArgs][nchar(arg[xtraArgs])>0]	
			#replace default args with defaults
			for(da in names(defArgs)){
				bdy = gsub(da, defArgs[da], bdy)
				#making the rest of the datastructure consistent
				if( is.symbol(defArgs[[da]]) ){
					eval(parse( text=sprintf("self$%s='%s'", da, defArgs[da]) ))
				}else{
					eval(parse( text=sprintf("self$%s=%s", da, defArgs[da]) ))
				}
			}
			#replace place base args with self variables
			for( ea in xtraArgs ){
		                #
		                var = sprintf("self$%s", ea)
		                bdy = gsub(ea, var, bdy)
		        }
			
		        #
		        eval(parse( text=sprintf("self$%s=function(%s){}", name, paste(mainArgs, collapse=',')) ))
		        eval(parse( text=sprintf("body(self$%s)=parse(text=c('{', bdy, '}'))", name) ))
		},
		#
                selfToPar = selfToPar,
                parToSelf = parToSelf
	)
)
