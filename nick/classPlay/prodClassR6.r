
#
library(R6)
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
		#mortality
		mn = NA,
		mf = NA,
		#functions
		AtoL = NA,
		LtoW = NA,
		#computational
		method = 'rk4',
		
		#
		initialize = function( 	N0   = NA,
					mn   = NA,
					mf   = NA,
					time = NA,  
					dNdt = NA	
		){
			#misc variable digestion
			self$N0 = N0
			self$mn = mn
			self$mf = mf
			self$time = time
			
			#dNdt
			stopifnot(is.function(dNdt))
			private$dNdt_classify(dNdt)	
			
			#preallocate N
			self$N0 = N0
			self$N  = matrix(NA, nrow=length(time), ncol=1)
			rownames(self$N) = sprintf("TIME %d", time)
		},
		#
		iterate = function(method=self$method){
			#prechecking 
			
			#digest and possibly change ode method	
			self$method = method	
		
			#preallocate N
                        self$N  = matrix(NA, nrow=length(self$time), ncol=1)
                        rownames(self$N) = sprintf("TIME %d", self$time)
			
			#solve 
        		self$N = ode(self$N0, self$time, private$dNdt, private$dNdt_par, method)[,2]
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




