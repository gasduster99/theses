
#
suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(GA, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(GGally, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(mvtnorm, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(deSolve, quietly=TRUE) ))



#
#CLASS
#

#
ageModel = R6Class("AgeModel", lock_objects=FALSE,
	#
	public = list(
		#pop
		N0 = NA,
		N  = NA,
		#time
		time = NA,
		TT   = NA,
		#max age
		A = NA,
		#model
		sdp = NA,
		sdo = NA,
		#normalization constant
                q = 1,
		model = list(observation="LN", process=NA),
                prior = list(),
		#functions
		SRR  = NA,
		AtoL = NA,
		LtoW = NA,
		#computational
		ODE_method = 'rk4',
		OPT_method = 'L-BFGS-B',
		
		#
		initialize = function( 	N0   = NA,	
					A    = NA,	
					time = NA,  
					rDev = NA,
					dNdt = NA, 
					SRR  = NA,
					AtoL = function(a,Linf,k,a0){ Linf*(1-exp(-k*(a-a0))) }, 
					LtoW = function(l,rho,psi){ rho*l^psi },
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
						
			#SRR
			stopifnot(is.function(SRR))
			private$classify(SRR)
			
			#AtoL
			stopifnot(is.function(AtoL))
			private$classify(AtoL)
			
			#LtoW
			stopifnot(is.function(LtoW))
			private$classify(LtoW)
			
			#preallocate N
			self$N0 = N0
			self$A  = A
			private$buildN(time=time, A=A)
			
			#recDevs
			self$rDev = rDev
				
		},
		
		##NOTE: this only functions for consectutive times starting at 1 to self$TT
		iterate = function(method=self$ODE_method){
			#prechecking 
			
			#digest ode method	
			self$ODE_method = method
			#digest new dNdt values
			private$dNdt_classify(dNdt)	
			#last minute allocation
			private$buildN(time=self$time, A=self$A)
			
			#
			Ws = self$LtoW(self$AtoL(self$As:self$A))
			
			#solve 
			timeStart = self$N0 + self$rDev[1]
        		self$N[1,] = ode(timeStart, self$time[1:self$A], private$dNdt, private$dNdt_par, method)[,2]
        		#copy the upper diagonal over
        		for(r in 2:self$A){ for(c in r:self$A){ self$N[r,c]=self$N[r-1,c] }}

        		#NOTE: this loop assumes times are consecutive starting at 1 to self$TT
        		i = 0
        		for(t in 1:(self$TT-2)){
        		        #Length of the diagonal which starts on row t+1 
        		        D = min(self$A,self$TT-t)
        		        timeStart = self$SRR(Ws%*%self$N[t,self$As:self$A]) + self$rDev[t+1]
				d = ode(timeStart, self$time[1:D], private$dNdt, private$dNdt_par, method)[,2]
				for(a in 1:D){
        		                self$N[t+a,a] = d[a]
        		                i = i+1
        		        }
        		}
        		self$N[self$TT,1] = self$SRR(Ws%*%self$N[t,self$As:self$A]) + self$rDev[self$TT]
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
                        #gaBoost  : a logical optionally (default F) indicating if a persistent 
                        #       genetic algorithm should be used to assist local optimization.
                        #       genetic algoithm repeates until first and second finite difference 
                        #       derivatives are successful and hessian is inverted. Optionally 
                        #       gaBoost may be given as a list containting names values for list(popSize, maxiter, run).
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
                        #NOTE: the likelihood only allows distributions parameterized in terms of its mean 
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
			#par = private$selfToPar(parNames)
                        #nome = names(par)
			#funny = function(x){ names(x)=nome; -fun(x) }
			#print(funny(par))
			
                        #possibly precondition guesses with ga
                        if( gaBoost[[1]]!=F ){
                                #
                                if(gaBoost[[1]]==T){
                                        gaBoost = list(
                                                popSize = 1e4,
                                                maxiter = 1e3,
                                                run     = 30
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
                                                        type    = 'real-valued',
                                                        fitness = function(x){ names(x)=nome; -fun(x) },
                                                        names   = nome,
                                                        lower   = lower,
                                                        upper   = upper,
                                                        popSize = gaBoost[['popSize']],
                                                        maxiter = gaBoost[['maxiter']],
                                                        run     = gaBoost[['run']],
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
                                                #self$plotMean()
                                                #self$plotBand()
                                                #points(data)
                                                print( err )
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
                        display = nome[!nome%in%c(".__enclos_env__", "initialize", "iterate", "optimize", "clone", "print", "model", "prior", "plotMean", "plotBand", "plotRS", "SRR", "LtoW", "AtoL")]
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
			parNome = c()
                        parValues = c()	
                        for(pn in parNames){
                                eval(parse( text=sprintf("parValues=c(parValues,self$%s)", pn) ))
                        	l = eval(parse( text=sprintf("length(self$%s)", pn) ))
				parNome = c(parNome, pn, rep("", l-1)) 
			}
			names(parValues) = parNome
                        #
                        return(parValues)
                },

                #NOTE: parValues should be passed with names
                parToSelf = function(parValues){
                        #check is names exist

                        #
                        parNames = names(parValues)
                        #
			i = 1
			
                        for(pn in parNames){
                                #update self
				if(pn!=""){
					nome = pn
					eval(parse( text=sprintf("self$%s=parValues[pn]", pn) ))
                        	} else{ 
					eval(parse( text=sprintf("self$%s=c(self$%s, parValues[i])", nome, nome) )) 
				}
				#
				i = i+1
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
	
		#NOTE: numPars may only be 1 presently
		classify = function(fun, numPars=1){
			#
			name = as.character(substitute(fun))
			arg = formals(fun)
			bdy = body(fun)[-1]
			extraArgs = names(arg[(numPars+1):length(arg)])
			for( ea in extraArgs ){
				#
				var = sprintf("self$%s_%s", name, ea)
				bdy = gsub(ea, var, bdy)
			}
			#NOTE: check the argument with numPars>1
			eval(parse( text=sprintf("self$%s=function(%s){}", name, names(arg[1:numPars])) ))
			eval(parse( text=sprintf("body(self$%s)=parse(text=c('{', bdy, '}'))", name) ))	
		},
		
		#
		buildN = function(time=NA, A=NA){ #, As=NA, Af=NA){
			#preallocate N
			if(!is.na(A) && !is.na(time)){
				#basic
				self$time = time
				self$TT = max(time)
				self$A  = A	
				self$N  = matrix(NA, nrow=length(time), ncol=A)
				colnames(self$N) = sprintf("AGE %d", 1:A)
				rownames(self$N) = sprintf("TIME %d", time)
				##fancy
				#if(!is.na(As) & !is.na(Af)){
				#	self$As = As
				#	self$Af = Af
        			#	colnames(self$N)[As] = sprintf("%s (As)", colnames(self$N)[As])
        			#	colnames(self$N)[Af] = sprintf("%s (Af)", colnames(self$N)[Af])
				#}
			}
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







##
##FUNCTIONS
##
#
##dNdt
#dNdt = function(t, y, par){
#        #t      : current time step as requested by ode
#        #y      : value at previous time as requested by ode
#        #par    : parameters to be passed by value
#                #af     : age suseptible to fishing given as integer
#                #mn     : natural mortality given as numeric
#                #mf     : fishing mortality given as numeric
#	
#        #unpack
#        af = par$af
#        mn = par$mn
#        mf = par$mf
#        #assume population is not in the fishery, but if it is include fishing mortality
#        out = y*(exp(-mn)-1)
#        if(t>=af){ out = y*(exp(-mn-mf)-1) }
#        #
#        return( list(out) )
#}
##SRR
#SRR = function(S, a, b, c){
#        #S : spawning biomass
#        #a : BH a parameter
#        #b : BH b parameter
#        #c : Shepard c parameter
#        #
#        #value : the number of recruits from the given biomass, assuming BH spawning
#        
#        #
#        return( a*S/(1+(b*S)^c) )
#}
#
##
##MAIN
##
#
##define functions
#am = ageModel$new( dNdt=dNdt, SRR=SRR, A=10, TT=200, Af=3, As=3 )
##SRR
#am$SRR_a = 3
#am$SRR_b = 1/50000
#am$SRR_c = 1
##AtoL
#am$AtoL_Linf = 50
#am$AtoL_a0   = 0
#am$AtoL_k    = 0.25
##LtoW
#am$LtoW_rho = 0.1
#am$LtoW_psi = 3
#
###define data structures
##am$A  = 10
##am$Af = 3
##am$As = 3
##am$TT = 200




