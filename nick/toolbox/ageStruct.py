import numpy as np
from scipy.integrate import solve_ivp

#
#

#
class ageModel(object):
	#
	def __init__(self, dNdt, SRR, AtoL, LtoW, cats):
		'''
		dNdt 	: mortality derivative given as a function
		SRR	: stock recruitment relationship given as a preloaded string or a function
		AtoL	: age to length relationship given as a function
		LtoW	: allometric length to weight relationship as a function 
		cats	: age/length categories given as an ordered numpy.array of ?strings?
		'''

		#
		self.dNdt = dNdt
		self.SRR  = SRR
		self.AtoL = AtoL
		self.LtoW = LtoW
		self.cats = cats
		#
		self.NN = np.array([])
	#
	
	#
	def iterate(self, N0, time, method):
		'''
		N0	: initial population size (virgin age zeros) given as a numeric
		time 	: numpy.array of requested times to interate through
		method	: ode iterator, given as a string, to be handed to scipy.integrate.ode
		'''
		
		#
		self.NN = np.empty([len(time),len(self.cats)], dtype=int)
		#
		t0 = min(time)
		tT = max(time)
		self.NN[0,:] = solve_ivp(dNdt, [t0,tT], N0, t_eval=time).y
		#odeObj.set_integrator(method)
		#odeObj.set_initial_values(N0)
		#odeObj.set_f_params()
		#odeObj.set_jac_params()
		
		#solve
		#N[1,] = ode(500000, time[1:A], dNdt, c(Af, nm, fm), method)[,2]
		
	#
	
	#
	
	
