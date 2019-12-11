rm(list=ls())

#
source('prodClassR6.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, mn, mf){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #mn     : natural mortality given as numeric
        #mf     : fishing mortality given as numeric
 
        #assume population is not in the fishery, but if it is include fishing mortality 
	out = y*(exp(-mn-mf)-1) 
        #
        return( list(out) )
}

#
#MAIN
#

#define functions
pm = prodModel$new( dNdt=dNdt, time=1:200, N0=500000, mn=0.2, mf=0.2 )
#AtoL
pm$AtoL_Linf = 50
pm$AtoL_a0   = 0
pm$AtoL_k    = 0.25
#LtoW
pm$LtoW_rho = 0.1
pm$LtoW_psi = 3

#
pm$iterate()


