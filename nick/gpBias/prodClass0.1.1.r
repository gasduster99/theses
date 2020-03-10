#
suppressWarnings(suppressMessages( library(R6, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(deSolve, quietly=TRUE) ))
#
source('modelFunk.r')
source('outputFunk.r')
source('classifyFunk.r')
source('parSelfParFunk.r')
source('optimizeFunk.r')

#
#CLASS
#

#
prodModel = R6Class("ProdModel", lock_objects=FALSE,
	#
	public = list(
		#pop
		N  = NA,
		N0 = NA,	
		N0Funk = NA,
		#time
		time = NA,
		#model
                lsdo = NA, 
		lq = 0, #log proportionality constant between cpue and N
		model = list(observation="LN"),
                prior = list(),
		#functions
		#computational
		ODE_method = 'rk4',
		OPT_method = 'L-BFGS-B',
		
		#
		initialize = function( 	N0   = NA,
					time = NA,
					dNdt = NA,
					N0Funk = NA,
					...
		){	
			#misc variable digestion
                        misc = list(...)
                        miscNames = names(misc)
                        for(m in miscNames){
                                eval(parse( text=sprintf("self$%s=misc[[m]]", m) ))
                        }
			
			#dNdt
			stopifnot(is.function(dNdt))
			private$dNdt_classify(dNdt)	
			
			#N0
			if( is.function(N0Funk) ){ 	
				private$N0_classify(N0Funk)
				N0 = self$N0Funk()
			}
			
			#preallocate N
			self$N0 = N0
			self$time = time
			self$N  = matrix(NA, nrow=length(time), ncol=1)
			rownames(self$N) = sprintf("TIME %d", time)
		},
		
		#
		iterate = function(method=self$ODE_method){
			#method	: optionally update the ODE method to be handed to the ode function (default 'rk4')
			#
			#value	: nothing returned, but self$N is updated with current values

			#prechecking 
			
			#digest and possibly change ode method	
			self$ODE_method = method	
	
			#N0
                        if( is.function(N0Funk) ){ N0=self$N0Funk() }
				
			#last minute allocation and variable updates
			self$N0 = N0
                        self$N = matrix(NA, nrow=length(self$time), ncol=1)
                        rownames(self$N) = sprintf("TIME %d", self$time)
			
			#solve 
        		self$N = ode(self$N0, self$time, private$dNdt, parms=NULL, method=method)[,2]
		},
		
		#
                optimize = optimize,
		
		#
		printer   = printSelf 
		printSelf = function(ins=c()){
			self$printer(ins, outs=c(
				"iterate", "optimize", "model",	"prior", 
				"plotMean", "plotBand", "plotRS", "N0Funk", 
				"save", "load", "printer"
			))
		},
		plotMean = plotMean,
		plotBand = plotBand,
		plotRS 	 = plotRS,
		#
		save = function(fileName){ saveRDS(self, file=fileName) },
		load = function(fileName){ readRDS(fileName) }
	),
	
	#
	private = list(
		##
		#dNdt = NA,
		
		#
		selfToPar = selfToPar,	
		parToSelf = parToSelf, #NOTE: parValues should be passed with names
		
		#
		dNdt_classify = dNdt_classify,
		N0_classify = N0_classify,
		classify = classify,
		
		#
		dLikes = dLikes,
		qLikes = qLikes
	)
)




