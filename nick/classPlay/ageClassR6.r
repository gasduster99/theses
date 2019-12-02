
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
		#
		dNdt = NA,
		SRR  = NA,
		AtoL = NA,
		LtoW = NA,
		#
		initialize = function( 	dNdt=NA, 
					SRR=NA,
					AtoL=function(a,Linf,k,a0){ Linf*(1-exp(-k*(a-a0))) }, 
					LtoW=function(l,rho,psi){ rho*l^psi }
		){
			#dNdt
			stopifnot(is.function(dNdt))
			self$dNdt = dNdt
			#SRR
			stopifnot(is.function(SRR))
			private$classify(SRR, 1)
			#AtoL
			stopifnot(is.function(AtoL))
			private$classify(AtoL, 1)
			#LtoW
			stopifnot(is.function(LtoW))
			private$classify(LtoW, 1)
		},
		iterate = function()
	),
	#
	private = list(
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
			#
			eval(parse( text=sprintf("self$%s=function(%s){}", name, names(arg[1:numPars])) ))
			eval(parse( text=sprintf("body(self$%s)=parse(text=bdy)", name) ))
			
			##
			#arg = formals(SRR)
			#bdy = body(SRR)[-1]
			#extraArgs = names(arg[2:length(arg)])
			#for( ea in extraArgs ){
			#	#
			#	var = sprintf("self$SRR_%s", ea)
			#	bdy = gsub(ea, var, bdy)
			#	#var = sprintf("SRR_%s", ea)
			#	#self$set("public", var, NA)
			#	#eval(parse( text=paste(var, "=NA") ))
			#}
			##
			#self$SRR = function(S){}
			#body(self$SRR) = parse(text = bdy)
			
		}
	)
)

#
#MAIN
#

#dNdt
dNdt = function(t, y, par){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #par    : parameters to be passed by value
                #af     : age suseptible to fishing given as integer
                #nm     : natural mortality given as numeric
                #fm     : fishing mortality given as numeric
        
        ##
        #par = list(...)
        #unpack
        af = par$af
        nm = par$nm
        fm = par$fm
        #assume population is not in the fishery, but if it is include fishing mortality
        out = y*(exp(-nm)-1)
        if(t>=af){ out = y*(exp(-nm-fm)-1) }
        #
        return( list(out) )
}
#SRR
SRR = function(S, a, b, c){
        #S : spawning biomass
        #...    : parameters to be passed by value
                #a : BH a parameter
                #b : BH b parameter
                #c : Shepard c parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

        ##
        #par = list(...)
        #a   = par$a
        #b   = par$b
        #c   = par$c
        
        #
        return( a*S/(1+(b*S)^c) )
}
a = 3
b = 1/50000
c = 1

#
am = ageModel$new( dNdt=dNdt, SRR=SRR )
#SRR
am$SRR_a = 3
am$SRR_b = 1/50000
am$SRR_c = 1
#AtoL
am$AtoL_Linf = 50
am$AtoL_a0   = 0
am$AtoL_k    = 0.25
#LtoW
am$LtoW_rho = 0.1
am$LtoW_psi = 0.1

