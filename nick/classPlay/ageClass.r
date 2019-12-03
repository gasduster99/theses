
#
library(methods)

#
#CLASS
#

#
ageModel = setClass("AgeModel", 
	slots = c(
		dNdt = "function",
		SRR  = "function",
		AtoL = "function",
		LtoW = "function",
		par  = "list"
	),
	prototype = list(
		AtoL = function(a){ par$Linf*(1-exp(-par$k*(a-par$a0))) },
		LtoW = function(l, rho){ rho*l^3 }
	)
)

#
#METHODS
#

#
#MAIN
#

#
am = ageModel()

#
am@SRR = function(S, ...){
        #S : spawning biomass
	#...    : parameters to be passed by value
        	#a : BH a parameter
        	#b : BH b parameter
		#c : Shepard c parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

	#
	par = list(...)
	a   = par$a
	b   = par$b
	c   = par$c
	
        #
        return( a*S/(1+(b*S)^c) )
}

#
am@dNdt = function(t, y, ...){
	#t	: current time step as requested by ode
	#y	: value at previous time as requested by ode
	#...	: parameters to be passed by value
		#af	: age suseptible to fishing given as integer
        	#nm 	: natural mortality given as numeric
        	#fm	: fishing mortality given as numeric
	
	#
	par = list(...)
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
am@par[['a']] = 3
am@par[['b']] = 1/500000
am@par[['c']] = 1
#dNdt
am@par[['af']] = 3
am@par[['nm']] = 0.2
am@par[['fm']] = 0.2
#AtoL
am@par[['Linf']] = 50
am@par[['a0']]   = 0
am@par[['k']]    = 0.25







                        ##
                        #arg = formals(SRR)
                        #bdy = body(SRR)[-1]
                        #extraArgs = names(arg[2:length(arg)])
                        #for( ea in extraArgs ){
                        #       #
                        #       var = sprintf("self$SRR_%s", ea)
                        #       bdy = gsub(ea, var, bdy)
                        #       #var = sprintf("SRR_%s", ea)
                        #       #self$set("public", var, NA)
                        #       #eval(parse( text=paste(var, "=NA") ))
                        #}
                        ##
                        #self$SRR = function(S){}
                        #body(self$SRR) = parse(text = bdy)

