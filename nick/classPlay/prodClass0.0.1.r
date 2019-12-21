
#
library(R6)
library(GA)
library(deSolve)

#
#CLASS
#

#
prodModel = R6Class("ProdModel", lock_objects=FALSE,
	#
	public = list(
		#pop
		N0 = NA,
		N  = NA,
		#time
		time = NA,
		#model
                sdp = NA,
                sdo = NA,
		#NOTE: the likelihood only distributions parameterized in terms of its mean and standard deviation
                #I need to reconsider how to build the likelihood/prior handling system
		q = 1,
		model = list(observation="LN", process=NA),
		#likelihood = list(observation=dlnorm, process=NA), 
                prior = list(), #parameter name=function
		#functions
		#computational
		ODE_method = 'rk4',
		OPT_method = 'L-BFGS-B',
		
		#
		initialize = function( 	N0   = NA,	
					time = NA,  
					dNdt = NA,
					...
		){	
			#misc variable digestion
                        misc = list(...)
                        miscNames = names(misc)
                        for(m in miscNames){
                                eval(parse( text=sprintf("self$%s=misc[[m]]", m) ))
                        }
			
			#dNdt
			stopifnot(is.function(dNdt))
			private$dNdt_classify(dNdt)	
	
			#preallocate N
			self$N0 = N0
			self$time = time
			self$N  = matrix(NA, nrow=length(time), ncol=1)
			rownames(self$N) = sprintf("TIME %d", time)
		},
		
		#
		iterate = function(method=self$ODE_method){
			#prechecking 
			
			#digest and possibly change ode method	
			self$method = method	
		
			#last minute allocation and variable updates
                        self$N  = matrix(NA, nrow=length(self$time), ncol=1)
                        rownames(self$N) = sprintf("TIME %d", self$time)
			
			#solve 
        		self$N = ode(self$N0, self$time, private$dNdt, parms=NULL, method=method)[,2]
       			#self$I = self$q*self$N
		},
		
		#
                optimize = function(    data,
					parNames,
                                        lower, upper,
                                        method=self$OPT_method,
                                        cov=F,
                                        gaBoost=F,
                                        control = list()
                ){
                        #prechecking

                        #digest opt method      
                        self$OPT_method = method	

                        #here I assemble the model
                        #NOTE: the likelihood only allows distributions parameterized in terms of its mean and standard deviation
                        fun = function(par){
				#unpack par
				private$parToSelf(par)
				#compute mean function	
				self$iterate()
				#evaluate likelihood
                                #like = self$likelihood$observation(data, self$q*self$N, self$sdo, log=T)
				like = private$likes[[self$model$observation]](self, data)
				##
				#print(self$q)
				#print(self$sdo)	
				#print(like)
				#print(self$N)
				##
				return( -sum(like) )
                        }
				
                        #possibly precondition guesses with ga
			if( gaBoost!=F ){
				#
				if(gaBoost==T){
					gaBoost = list(
						popSize = 1e4,
						maxiter = 1e3,
						run	= 30
					)
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
						out = list()
						#
						out[['gaOut']] = ga(
							type 	= 'real-valued', 
							fitness	= function(x){ names(x)=nome; -fun(x) },
							names	= nome,
							lower	= lower,
						        upper   = upper,
						       	popSize = gaBoost[['popSize']], 
						       	maxiter = gaBoost[['maxiter']],
							run	= gaBoost[['run']],
						       	optim   = T,
							parallel= T,
							#monitor = F,
							suggestions = par
						)
						#make sure to update with the best solution
						private$parToSelf(out[['gaOut']]@solution[1,])
						par = out[['gaOut']]@solution
	
						#optim
                        			out[['optimOut']] = optim(
                        			        par[1,],
                        			        fun,
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

                        			#?how to handle covariance?
                        			if( cov ){ self$rsCov = solve(out$optimOut$hessian) }
						
						#
						flag = F
					}, error=function(err){
						#
						writeLines( sprintf("\nRound: %s\n", i) ) 
						print( self )
						#print( err )
						#
						return(T) 
					}
					)
					#
					i = i+1
				}

			} else{
                        	#optim
                        	out = optim(
					private$selfToPar(parNames), 
					fun,
                        		lower = lower,
                        		upper = upper,
                        		hessian = cov,
					method  = self$OPT_method, 
                        		control = control
                        	)	

				#?how to handle covariance?
				if( cov ){ self$rsCov = solve(out$optimOut$hessian) }
			}
			#
			return( out )
		},
		
		#
		print = function(){
			#
			nome = names(self)
			display = nome[!nome%in%c(".__enclos_env__", "initialize", "iterate", "optimize", "clone", "print")]
			#
			for(d in display){
				#
				writeLines(sprintf('%s:', d))
				typeof( eval(parse( text=sprintf("self$%s", d) )) )
			}
		}
		
		##
		#plot = function(){
		#}	
	),
	
	#
	private = list(
		#
		dNdt = NA,
		
		#
		selfToPar = function(parNames){
			#check if variable names exist

			#
			parValues = c()
			for(pn in parNames){
				eval(parse( text=sprintf("parValues['%s']=self$%s", pn, pn) ))
			}
			#
			return(parValues)
		},
		
		#NOTE: parValues should be passed with names
		parToSelf = function(parValues){
			#check is names exist

			#
			parNames = names(parValues)
			#
			for(pn in parNames){
				#update self
				eval(parse( text=sprintf("self$%s=parValues[pn]", pn) ))	
			}
		},
		
		#
		dNdt_classify = function(fun, numPars=2){
			#
                        name = as.character(substitute(fun))
                        arg = formals(fun)
                        bdy = body(fun)[-1]
                        mainArgs = names(arg[1:(numPars)])
                        xtraArgs = names(arg[(numPars+1):length(arg)])
                        for( ea in xtraArgs ){
                                #NOTE: i need to identify variable better from coincidental text
                                var = sprintf("self$%s", ea)
                                bdy = gsub(ea, var, bdy)
                        }
			#
			eval(parse( text=sprintf("private$dNdt=function(%s, %s, dNdt_par){}", names(arg[1]), names(arg[2])) ))
                        body(private$dNdt) = parse( text=c('{', bdy, '}') )
		},
		
		#
		classify = function(fun, numPars=1){
			#
			
			#
			name = as.character(substitute(fun))
			arg = formals(fun)
			bdy = body(fun)[-1]
			mainArgs = names(arg[1:(numPars)])
			xtraArgs = names(arg[(numPars+1):length(arg)])
			for( ea in xtraArgs ){
				#
				var = sprintf("self$%s_%s", name, ea)
				bdy = gsub(ea, var, bdy)
			}
			
			#
			eval(parse( text=sprintf("self$%s=function(%s){}", name, paste(mainArgs, collapse=',')) ))
			eval(parse( text=sprintf("body(self$%s)=parse(text=c('{', bdy, '}'))", name) ))	
		},
		
		#
		likes = list(
			#
			LN = function(self, data){
				dnorm(log(data), log(self$q)+log(self$N), self$sdo, log=T)
			},
			#
			N = function(self, data){
				dnorm(data, self$q*self$N, self$sdo, log=T)
			}
		)
	)
)




