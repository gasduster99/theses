
#
library(R6)
library(deSolve)

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
		likelihood = list(observation=dlnorm, process=NA), #NOTE: the likelihood only takes normal and lognormal observation and process models
		prior = list(),#parameter name=function
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
		optimize = function(	pars, 
					lower, upper, 
					method=self$OPT_method, 
					cov=F, 
					gaBoost=F, 
					control = list()
		){
			#pars is a vector of variable names 
			#estimates update self variables, and an optional 
			#covariance matrix is defined (or returned?)
			
			#prechecking
			
			#digest opt method      
                        self$OPT_method = method
			
			#I need to build the likelihood/prior handling system
			
			#build par vector from self
			
			#here I assemble the model
			#NOTE: the likelihood only takes normal and lognormal observation and process models
			fun = function(){
				self$likelihood$observation()
			}
			
			#possibly precondition guesses with ga
			
			#optim
			optimOut = optim(par, fun, 
				lower = lower, 
				upper = upper,
				hessian = cov,
				control = control
			)
		}
	),
	#
	private = list(
		#
		dNdt = NA,
		dNdt_par = c(),	
		#
		dNdt_classify = function(fun){
			#
			arg = formals(fun)
			bdy = as.character(body(fun)[-1])
			extraArgs = names(arg[3:length(arg)])
			parBlock = c()
			for( ea in extraArgs ){
				#fill dNdt_par
				eval(parse( text=sprintf("private$dNdt_par[ea]=self$%s", ea) ))
				#build par block
				parBlock = c(parBlock, sprintf("%s = dNdt_par['%s']", ea, ea))
			}
			#
			eval(parse( text=sprintf("private$dNdt=function(%s, %s, dNdt_par){}", names(arg[1]), names(arg[2])) ))
			body(private$dNdt) = parse( text=c('{', parBlock, bdy, '}') )
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
		}
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




