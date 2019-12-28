#
suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(GA, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(GGally, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(mvtnorm, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(deSolve, quietly=TRUE) ))

#
#FUNCTIONS
#

#
makeTransparent = function(someColor, alpha=100){
	newColor = col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){
		rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3], alpha=alpha, maxColorValue=255)
	})
}

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
                sdp = NA, #NOTE: Do I want this here, or possibly only when I create a model with process uncertainty
                sdo = NA,
                #normalization constant
		q = 1,
		model = list(observation="LN", process=NA),
                prior = list(),
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
			#method	: optionally update the ODE method to be handed to the ode function (default 'rk4')
			#
			#value	: nothing returned, but self$N is updated with current values

			#prechecking 
			
			#digest and possibly change ode method	
			self$ODE_method = method	
		
			#last minute allocation and variable updates
                        self$N = matrix(NA, nrow=length(self$time), ncol=1)
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
			#data	  : a vector of data used to fit specified model.
			#	In the future consider adding data to object as an inheritted form.
			#parNames : a vector of strings matching the names of the parameters to be 
			#	optimized.
			#lower	  : a vector of lower bounds for searching parameter space
			#upper    : a vector of upper bounds for searching parameter space
			#method   : a string optionally defining the method of optimization to be 
			#	handed to optim (default 'L-BFGS-B')
			#cov	  : a logical optionally indicating if hessian should be computed and 
			#	inverted in optimization process
			#gaBoost  : a logical optionally (default F) indicating if a persistent 
			#	genetic algorithm should be used to assist local optimization.
			#	genetic algoithm repeates until first and second finite difference 
			#	derivatives are successful and hessian is inverted. Optionally 
			#	gaBoost may be given as a list containting names values for list(popSize, maxiter, run).
			#control  : additional control parameters to be passed to optim
			#
			#value    : optimization objects are returned. Parameters values are 
			#	updated inplace. rsCov is defined to self.
			
                        #prechecking

                        #digest opt method      
                        self$OPT_method = method
			
			#
			out = list()

                        #here I assemble the model
                        #NOTE: the likelihood only allows distributions parameterized in terms of its mean and standard deviation
                        fun = function(par){
				#unpack par
				private$parToSelf(par)
				#compute mean function	
				self$iterate()
				#evaluate likelihood
				like = private$dLikes[[self$model$observation]](self, data)
				#	
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
                        	out[['optimOut']] = optim(
					private$selfToPar(parNames), 
					fun,
                        		lower = lower,
                        		upper = upper,
                        		hessian = cov,
					method  = self$OPT_method, 
                        		control = control
                        	)

				#how to handle covariance
				if( cov ){ self$rsCov = solve(out$optimOut$hessian) }
			}
			#
			return( out )
		},
		
		#
		print = function(){
			#
			n = 5
			#
			nome = names(self)
			display = nome[!nome%in%c(".__enclos_env__", "initialize", "iterate", "optimize", "clone", "print", "model", "prior", "plotMean", "plotBand", "plotRS")]
			display = display[order(nchar(display))]
			#
			for(d in display){
				#
				text = sprintf("self$%s", d)
				if( typeof(eval(parse(text=text)))=='list' ){
					#
					cat( sprintf('%s\t:\n', d) )
					print( eval(parse(text=text)) )
				} else{
					#	
					if(length(eval(parse( text=text )))>n){
						cat( sprintf('%s\t:', d), eval(parse( text=text ))[1:n], '...\n' )
					} else{
						cat( sprintf('%s\t:', d), eval(parse( text=text )), '\n' )
					}
				}
			}
		},
		
		#
		plotMean = function(col='black', alpha=100, lwd=3, add=F){
			#arguments passed directly to plotting functions
			
			#	
                        if(!add){
                                plot(self$time, self$q*self$N, 
                                        type='l',
                                        lwd=lwd, 
                                        col=col
                                )
                        } else{
                                lines(self$time, self$q*self$N,
                                        lwd=lwd,
                                        col=col
                                )
                        }

		},
		
		#
		plotBand = function(prob=0.95, col='black', alpha=100){
			#arguments passed directly to plotting functions
			
			#
			left = (1-prob)/2
			right = prob+left
			#
			polygon( c(self$time, rev(self$time)), 
				c(
					private$qLikes[[self$model$observation]](self, left), 
					rev(private$qLikes[[self$model$observation]](self, right))
				), 
				col = makeTransparent(col, alpha=alpha), 
				border = NA
			)	
		},
		
		#
		plotRS = function(m=10^4, save=F){
			#m	: how many samples
			#save	: FALSE or a filename
			
			#
			parNames = colnames(self$rsCov)	
			sam = rmvnorm(m, private$selfToPar(parNames), self$rsCov)
			#
			if(save==F){
				 GGally::print_if_interactive(ggpairs(as.data.frame(sam)))
			} else{	
				ggsave(filename=save, plot=ggpairs(as.data.frame(sam)))
			}
			#
			return(sam)
		}	
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
		dLikes = list(
			#
			LN = function(self, data){
				dnorm(log(data), log(self$q)+log(self$N), self$sdo, log=T)
			},
			#
			N = function(self, data){
				dnorm(data, self$q*self$N, self$sdo, log=T)
			}
		),
		
		#
		qLikes = list(
			#
			LN = function(self, prob){
				qlnorm(prob, log(self$q)+log(self$N), self$sdo)
			},
			#
			N = function(self, prob){
				qnorm(prob, self$q*self$N, self$sdo)
			}
		)
	)
)




