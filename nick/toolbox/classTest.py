#!/usr/bin/env python3

#
import numpy as np
#
from ageStruct import *

#
#FUNCTIONS
#

#must start with t, y, followed by *args or **kwargs
def dNdt(t, y, Af, nm, fm):
	#
	out = y*(exp(-nm)-1)
	if t>=Af: out = y*(exp(-nm-fm)-1)
	#
	return( out )
#

#
cats = np.arange(10)
parms = {
	'Af' : 3,
	'nm' : 0.2,
	'fm' : 0.2
}
mod  = ageModel(dNdt, 0, 0, 0, cats, parms)

