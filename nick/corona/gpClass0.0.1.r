#
suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(GA, quietly=TRUE) ))
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
		OPT_method = 'L-BFGS-B',
		#
		initialize = function(
			X  = NA,
			Y  = NA,
			B  = NA,
			S2 = NA,
			...
		){
			#misc variable digestion
                        misc = list(...)
                        miscNames = names(misc)
                        for(m in miscNames){
                                eval(parse( text=sprintf("self$%s=misc[[m]]", m) ))
                        }
			
			##
			#self$yRes = self$Y-predict(self$lm)
			
			#some check of B, X, and Y	
			self$X = X
			self$Y = Y
			self$B = B

			#
			if(!'obsV'%in%miscNames){ self$obsV=rep(0, length(self$Y)) }
			
			#Covariance function
                        stopifnot(is.function(S2))
                        private$classify(S2, 2)
		},
			
		#
		fit = function(	parNames,
                        	lower, upper,
				method=self$OPT_method,
                        	cov=F,
				gaBoost=F,
                        	control=list()
		){
			#digest opt method      
       			self$OPT_method = method

       		 	#
       		 	out = list()
			
			#
			loglikeL = function(par){
			        #
				private$parToSelf(par)
				#
				X = cbind(1, self$X)
				SIG2 = self$S2(self$X, self$X) + diag(self$obsV) + self$g*diag(length(self$obsV))
				#-dmvnorm(self$yRes, sigma=SIG2, log=T)
				-dmvnorm(self$Y, X%*%self$B, sigma=SIG2, log=T)
			}

	        	#possibly precondition guesses with ga
	        	if( gaBoost[[1]]!=F ){
	        	        #
	        	        set = list(
	        	                popSize = 1e4,
	        	                maxiter = 1e3,
	        	                run     = 30,
	        	                parallel= T
	        	        )
	        	        for(bn in names(gaBoost)){
	        	                set[[bn]] = gaBoost[[bn]]
	        	        }
				
				#
				par = private$selfToPar(parNames)
                		nome = names(par)
				
                		# 
                		i = 0
                		flag = T
		                while( flag ){
		                        tryCatch({
		                                #
		                                out[['gaOut']] = ga(
		                                        type    = 'real-valued',
		                                        fitness = function(x){ names(x)=nome; -loglikeL(x) },
		                                        names   = nome,
		                                        lower   = lower,
		                                        upper   = upper,
		                                        popSize = set[['popSize']], #gaBoost[['popSize']], 
		                                        maxiter = set[['maxiter']], #gaBoost[['maxiter']],
		                                        run     = set[['run']],     #gaBoost[['run']],
		                                        optim   = T,
		                                        parallel= set[['parallel']],        #T,
		                                        #monitor = F, 
		                                        suggestions = par
		                                )
		                                #make sure to update with the best solution
		                                private$parToSelf(out[['gaOut']]@solution[1,])
		                                par = out[['gaOut']]@solution
		
		                                #optim
		                                out[['optimOut']] = optim(
		                                        par[1,],
		                                        loglikeL,
		                                        lower = lower,
		                                        upper = upper,
		                                        hessian = cov,
		                                        method  = self$OPT_method,
		                                        control = control
		                                )
		
		                                #make sure to update with the best solution
		                                if( out[['gaOut']]@fitnessValue>out[['optimOut']]$value ){
		                                        #
		                                        private$parToSelf(out[['gaOut']]@solution[1,])
		                                }
		
		                                #
		                                if( cov ){
		                                        #
		                                        self$rsCov = chol2inv(chol(out$optimOut$hessian))
		                                        colnames(self$rsCov) = nome
		                                        rownames(self$rsCov) = nome
		                                }
		
		                                #
		                                flag = F
		                        }, error=function(err){	
		                                #
		                                writeLines( sprintf("\nRound: %s\n", i) )
		                                self$printSelf()
		                                print(err) 
		                                #
		                                return(T)
		                        }
		                        )
		                        # 
		                        i = i+1
		                }
        		} else{
				#
				par = private$selfToPar(parNames)
                                nome = names(par)
				
				#
				out[['optimOut']] = optim(
                        		par,
                        		loglikeL,
                        		lower = lower,
                        		upper = upper,
                        		hessian = cov,
                        		method  = 'L-BFGS-B',
                        		control = control
                		)
				#
				self$aic = 2*(out[['optimOut']]$value + length(parNames))
				#
                		if( cov ){
                		        #
                		        self$rsCov = chol2inv(chol(out$optimOut$hessian))
                		        colnames(self$rsCov) = nome
                		        rownames(self$rsCov) = nome
                		}
			}
			#
			return( out )
		},
		
		#
		predictMean = function(XX){
			if( !'KInv' %in% ls(self) ){	
				self$KInv = chol2inv(chol( self$S2(self$X, self$X) + diag(self$obsV) + self$g*diag(length(self$obsV)) ))
			}
			Xint  = cbind(1, self$X) #as.matrix(cbind(1, self$X))
			XXint = cbind(1, XX)	 #as.matrix(cbind(1, XX))
			#predict(self$lm, XX) + self$S2(self$X, XX)%*%self$KInv%*%self$yRes
			XXint%*%self$B + self$S2(self$X, XX)%*%self$KInv%*%(self$Y-Xint%*%self$B)
		},

		#
		predictVar = function(XX){
			if( !'KInv' %in% ls(self) ){	
				self$KInv = chol2inv(chol( self$S2(self$X, self$X) + diag(self$obsV) + self$g*diag(length(self$obsV)) ))
			}
			#KPP-KPX%*%KInv%*%KXP
			self$S2(XX, XX)	- self$S2(self$X, XX)%*%self$KInv%*%self$S2(XX, self$X)
		},
	
		#
		printer = printSelf,
                printSelf = function(ins=c()){
                        self$printer(ins, outs=c(
                                "S2", "fit", "save", "load", "lm", "printer", 
				"predictMean", "predictVar"
                        ))
                },
		
		#
                save = function(fileName){ saveRDS(self, file=fileName) },
                load = function(fileName){ readRDS(fileName) }
	),
	#
	private = list(
		#NOTE: maybe update classifyFunk with this function
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
