rm(list=ls())

#
source('ageClassR6.r')

#
#FUNCTIONS
#

#dNdt
dNdt = function(t, y, par){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #par    : parameters to be passed by value
                #af     : age suseptible to fishing given as integer
                #mn     : natural mortality given as numeric
                #mf     : fishing mortality given as numeric

        #unpack
        Af = par['Af']
        mn = par['mn']
        mf = par['mf']
        #assume population is not in the fishery, but if it is include fishing mortality
        out = y*(exp(-mn)-1)
        if(t>=Af){ out = y*(exp(-mn-mf)-1) }
        #
        return( list(out) )
}
dNdt = function(t, y, Af, mn, mf){
        #t      : current time step as requested by ode
        #y      : value at previous time as requested by ode
        #Af     : age suseptible to fishing given as integer
        #mn     : natural mortality given as numeric
        #mf     : fishing mortality given as numeric
 
        #assume population is not in the fishery, but if it is include fishing mortality
        out = y*(exp(-mn)-1)
        if(t>=Af){ out = y*(exp(-mn-mf)-1) }
        #
        return( list(out) )
}

#SRR
SRR = function(S, a, b, c){
        #S : spawning biomass
        #a : BH a parameter
        #b : BH b parameter
        #c : Shepard c parameter
        #
        #value : the number of recruits from the given biomass, assuming BH spawning

        #
        return( a*S/(1+(b*S)^c) )
}

#
#MAIN
#

#define functions
am = ageModel$new( dNdt=dNdt, SRR=SRR, time=1:200, N0=500000, A=10, Af=3, As=3, mn=0.2, mf=0.2 )
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
am$LtoW_psi = 3

#N0 = 500000
##define data structures
#am$A  = 10
#am$Af = 3
#am$As = 3
#am$TT = 200

##
#am$iterate()


