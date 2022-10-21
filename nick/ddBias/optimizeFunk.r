#
suppressWarnings(suppressMessages( library(GA, quietly=TRUE) ))

#
#

#
optimize = function(    data,
                        parNames,
                        lower, upper,
                        method=self$OPT_method,
                        cov=F,
			fitQ=T,
                        gaBoost=F,	
			persistFor=50,
                        control = list()
){      
        #data     : a vector of data used to fit specified model.
        #       In the future consider adding data to object as an inheritted form.
        #parNames : a vector of strings matching the names of the parameters to be 
        #       optimized.
        #lower    : a vector of lower bounds for searching parameter space
        #upper    : a vector of upper bounds for searching parameter space
        #method   : a string optionally defining the method of optimization to be 
        #       handed to optim (default 'L-BFGS-B')
        #cov      : a logical optionally indicating if hessian should be computed and 
        #       inverted in optimization process
	#fitQ	  : a logical indicating whether to find the MLE of log(q):lq via profile MLE,
	#	or if False don't change the initialized value of lq.
        #gaBoost  : a logical optionally (default F) indicating if a persistent 
        #       genetic algorithm should be used to assist local optimization.
        #       genetic algoithm repeates until first and second finite difference 
        #       derivatives are successful and hessian is inverted. Optionally 
        #       gaBoost may be given as a list containting names values for list(popSize, maxiter, run).
        #persistFor : a numeric indicating how many iterations of optimization tryCatch to engange in.
	#	set to Inf
	#control  : additional control parameters to be passed to optim
        #
        #value    : optimization objects are returned. Parameters values are 
        #       updated inplace. rsCov is defined to self.
        
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
		#check for pathological series
		if( any(is.na(self$N)) ){ return(NA) }
		if( any(is.na(self$B)) ){ return(NA) }
		if( any(!self$N>0) ){ return(NA) }
		if( any(!self$B>0) ){ return(NA) }
		#if( self$N[1]==0 ){ return(NA) }
		#profile maximization of log(q)
		if( fitQ ){ self$lq=mean(log(data)-log(self$B)) }
                #evaluate likelihood
                like = private$dLikes[[self$model$observation]](self, data) 
                return( -sum(like) )
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
                #if(gaBoost[[1]]==T){
                #        gaBoost = list( 
                #                popSize = 1e4,
                #                maxiter = 1e3,
                #                run     = 30
                #        )
                #}
                #
                par = private$selfToPar(parNames)
                nome = names(par)	
		
                # 
                i = 0
                flag = T       
                while( flag ){ 
                        flag = tryCatch({
                                #
                                out[['gaOut']] = ga(
                                        type    = 'real-valued', 
                                        fitness = function(x){ names(x)=nome; -fun(x) },
                                        names   = nome,
                                        lower   = lower,
                                        upper   = upper,
                                        popSize = set[['popSize']], #gaBoost[['popSize']], 
                                        maxiter = set[['maxiter']], #gaBoost[['maxiter']],
                                        run     = set[['run']],     #gaBoost[['run']],
                                        optim   = T,
                                        parallel= set[['parallel']], 	    #T,
                                        #monitor = F, 
                                        suggestions = par
                                )
                                #make sure to update with the best solution
                                private$parToSelf(out[['gaOut']]@solution[1,])
                                par = out[['gaOut']]@solution
				#profile maximization of log(q)
				if( fitQ ){
					self$iterate()
					self$lq=mean(log(data)-log(self$N)) 
				}
				                
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
					#profile maximization of log(q)
					if( fitQ ){ 
						self$iterate()
						self$lq=mean(log(data)-log(self$N)) 
					}
                                }
                                
                                #?how to handle covariance?
                                if( cov ){ 
					#
					self$rsCov = chol2inv(chol(out$optimOut$hessian)) 
					colnames(self$rsCov) = parNames
					rownames(self$rsCov) = parNames
				}
                                
                                #flag = F
				return( F )
                        }, error=function(err){
                                #
                                writeLines( sprintf("\nRound: %s\n", i) )
                                self$printSelf()
				print(err)
				
				#if(!is.na(self$N[1])){
                                #	self$plotMean()
                                #	self$plotBand()
				#}
                                #points(data)
                                #do.call(cat, err)
                                
				#be persistent but potentially not infinitely soo
                                #if(i<=persistFor){ return(T) }else{ return(F) }
                        	#flag = i<=persistFor
				return(i<=persistFor)
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
		if( cov ){ 
			#
			self$rsCov = chol2inv(chol(out$optimOut$hessian)) 
			colnames(self$rsCov) = parNames
			rownames(self$rsCov) = parNames
		}
                #if( cov ){ self$rsCov = solve(out$optimOut$hessian) }
        }
        #
        return( out )
}
