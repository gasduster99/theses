rm(list=ls())

library(methods)

#
#
#

#
ball = setClass('ball', 
	slots = c(
		color = "character",
		size = "numeric"
	)
)

setGeneric("deflate", function(b){standardGeneric("deflate")})
setMethod("deflate", 'ball',
	function(b){ 
		b@size = b@size-1 
		return( b )	
	}
)

b1 = ball(color='red', size=10)
show(b1)
b1 = deflate(b1)
show(b1)
