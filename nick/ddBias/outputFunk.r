#
suppressWarnings(suppressMessages( library(GGally, quietly=TRUE) ))
suppressWarnings(suppressMessages( library(mvtnorm, quietly=TRUE) ))

#
#

#
makeTransparent = function(someColor, alpha=100){
        newColor = col2rgb(someColor)
        apply(newColor, 2, function(curcoldata){
                rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3], alpha=alpha, maxColorValue=255)
        })
}

#
#

#
printSelf = function(ins, outs=c()){
        #
        n = 5
        #
        nome = names(self)
	extraNames = c(outs, ".__enclos_env__", "initialize", "printSelf", "clone")#, "iterate", "optimize", "clone", "model", "prior", "plotMean", "plotBand", "plotRS", "N0Funk", "B0Funk", "save", "load", "like")
        display = nome[!nome%in%extraNames | nome%in%ins]
	display = display[order(nchar(display))]
	if(length(ins)>0){ display=display[display%in%ins] } 
        #
        for(d in display){
                #
                text = sprintf("self$%s", d)
                if( typeof(eval(parse(text=text)))=='list' ){
                        #
                        cat( sprintf('%s\t:\n', d) )
                        print( eval(parse(text=text)) )
                } else{
                        #       
                        if(length(eval(parse( text=text )))>n){
                                cat( sprintf('%s\t:', d), eval(parse( text=text ))[1:n], '...\n' )
                        } else{
                                cat( sprintf('%s\t:', d), eval(parse( text=text )), '\n' )
                        }
                }
        }
}

#
plotQuan = function(quan=function(B){B}, col='black', alpha=100, lwd=3, add=F, ...){
        #arguments passed directly to plotting functions
	
	#
	private$N0_classify(quan)

        #       
        if(!add){
                plot(self$time, self$quan(),
                        type='l',
                        lwd=lwd,
                        col=col,
			...
                )
        } else{
                lines(self$time, self$quan(),
                        lwd=lwd,
                        col=col,
			...
                )
        }
}

#
plotMean = function(col='black', alpha=100, lwd=3, add=F){
        #arguments passed directly to plotting functions

        #       
        if(!add){
                plot(self$time, exp(self$lq+log(self$B)),
                        type='l',
                        lwd=lwd,
                        col=col
                )
        } else{
                lines(self$time, exp(self$lq+log(self$B)),
                        lwd=lwd,
                        col=col
                )
        }
}

#
plotBand = function(prob=0.95, col='black', alpha=100){
        #arguments passed directly to plotting functions

        #
        left = (1-prob)/2
        right = prob+left
        #
        polygon( c(self$time, rev(self$time)),
                c(
                        private$qLikes[[self$model$observation]](self, left),
                        rev(private$qLikes[[self$model$observation]](self, right))
                ),
                col = makeTransparent(col, alpha=alpha),
                border = NA
        )
}

#
plotRS = function(m=10^4, sample=F, save=F){
        #m      : how many samples
        #save   : FALSE or a filename

        #
        parNames = colnames(self$rsCov)
        sam = rmvnorm(m, private$selfToPar(parNames), self$rsCov)
        #
        if(save==F){
                GGally::print_if_interactive(ggpairs(as.data.frame(sam)))
        } else{
                ggsave(filename=save, plot=ggpairs(as.data.frame(sam)))
        }
        #
        if(sample){ return(sam) }
}
