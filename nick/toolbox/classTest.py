#!/usr/bin/env python3

#
import numpy as np
#
from ageStruct import *

#
#FUNCTIONS
#

#must start with t, y, followed by *args or **kwargs
def dNdt(t, y, af, nm, fm):
	#
	out = y*(exp(-nm)-1)
        if t>=af: out = y*(exp(-nm-fm)-1)
	#
	return( out )
#

##
#cats = np.arange(10)
#mod  = ageModel(0, 0, 0, 0, cats)

