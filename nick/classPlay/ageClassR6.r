
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
		#age
		A  = NA,
		As = NA,
		Af = NA,	
		#mortality
		mn = NA,
		mf = NA,
		#functions
		dNdt = NA,
		SRR  = NA,
		AtoL = NA,
		LtoW = NA,
		#computational
		method = 'rk4',
		#
		initialize = function( 	N0   = NA,	
					A    = NA,	
					As   = NA,
					Af   = NA,
					mn   = NA,
					mf   = NA,
					time = NA,  
					dNdt = NA, 
					SRR  = NA,
					AtoL = function(a,Linf,k,a0){ Linf*(1-exp(-k*(a-a0))) }, 
					LtoW = function(l,rho,psi){ rho*l^psi }
		){
			#misc variable digestion
			self$N0 = N0
			self$A  = A
			self$Af = Af
			self$As = As
			self$mn = mn
			self$mf = mf
			
			#dNdt
			stopifnot(is.function(dNdt))
			self$dNdt = dNdt
			private$dNdt_classify(dNdt)
			
			
			#SRR
			stopifnot(is.function(SRR))
			private$classify(SRR, 1)
			
			#AtoL
			stopifnot(is.function(AtoL))
			private$classify(AtoL, 1)
			
			#LtoW
			stopifnot(is.function(LtoW))
			private$classify(LtoW, 1)
			
			#preallocate N
			self$N0 = N0
			private$buildN(time=time, A=A, As=As, Af=Af)
			
			
		},
		#
		iterate = function(method=self$method){
			#prechecking 
			
			#digest and possibly change ode method	
			self$method = method	
			
			#last minute allocation
			private$buildN(time=self$time, A=self$A, As=self$As, Af=self$Af)
			
			Ws = self$LtoW(self$AtoL(self$As:self$A))
			
			#solve
        		self$N[1,] = ode(self$N0, self$time[1:self$A], self$dNdt, c(Af=self$Af, mn=self$mn, mf=self$mf), method)[,2]
        		#copy the upper diagonal over
        		for(r in 2:self$A){ for(c in r:self$A){ self$N[r,c]=self$N[r-1,c] }}

        		#NOTE: this loop assumes times are consecutive starting at 1 to self$TT
        		i = 0
        		for(t in 1:(self$TT-2)){
        		        #Length of the diagonal which starts on row t+1 
        		        D = min(self$A,self$TT-t);
        		        d = ode(self$SRR(Ws%*%self$N[t,self$As:self$A]), self$time[1:D], self$dNdt, c(Af=self$Af, mn=self$mn, mf=self$mf), method)[,2]
        		        for(a in 1:D){
        		                self$N[t+a,a] = d[a]
        		                i = i+1
        		        }
        		}
        		self$N[self$TT,1] = self$SRR(Ws%*%self$N[t,self$As:self$A]) 
		}
	),
	#
	private = list(
		#
		dNdt_par = list(),
		#
		dNdt_classify = function(fun){
			#
			arg = formals(fun)
			bdy = body(fun)[-1]
			extraArgs = names(arg[3:length(arg)])
			for( ea in extraArgs ){
				#fill dNdt_par
				eval(parse( text=sprintf("private$dNdt_par[ea]=self$%s", ea) ))
				#build par block
			}
		},
		#	
		classify = function(fun, numPars){
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
			print(bdy)
			#
			eval(parse( text=sprintf("self$%s=function(%s){}", name, names(arg[1:numPars])) ))
			eval(parse( text=sprintf("body(self$%s)=parse(text=bdy)", name) ))	
		},
		#
		buildN = function(time=NA, A=NA, As=NA, Af=NA){
			#preallocate N
			if(!is.na(A) && !is.na(time)){
				#basic
				self$time = time
				self$TT = max(time)
				self$A  = A	
				self$N  = matrix(NA, nrow=length(time), ncol=A)
				colnames(self$N) = sprintf("AGE %d", 1:A)
				rownames(self$N) = sprintf("TIME %d", time)
				#fancy
				if(!is.na(As) & !is.na(Af)){
					self$As = As
					self$Af = Af
        				colnames(self$N)[As] = sprintf("%s (As)", colnames(self$N)[As])
        				colnames(self$N)[Af] = sprintf("%s (Af)", colnames(self$N)[Af])
				}
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




