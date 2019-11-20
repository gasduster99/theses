#!/usr/bin/env python3

#
import numpy as np
#
from ageStruct import *

#
#FUNCTIONS
#

#
def dNdt(t, y, params):
	#
	
#


dNdt = function(t, y, params){
        #params[1] : af age suseptible to fishing
        #params[2] : nm natural mortality
        #params[3] : fm fishing mortality

        #unpack
        af = params[1]
        nm = params[2]
        fm = params[3]
        #assume population is not in the fishery, but if it is include fishing mortality
        out = y*(exp(-nm)-1)
        if(t>=af){ out = y*(exp(-nm-fm)-1) }
        #
        return( list(out) )
}



#
cats = np.arange(10)
mod  = ageModel(0, 0, 0, 0, cats)

