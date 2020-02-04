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
printSelf = function(){
        #
        n = 5
        #
        nome = names(self)
	extraNames = c(
	".__enclos_env__", "initialize", "iterate", "optimize", "clone", "printSelf", 
	"model", "prior", "plotMean", "plotBand", "plotRS", "N0Funk"
	)
        display = nome[!nome%in%extraNames]
        display = display[order(nchar(display))]
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
plotMean = function(col='black', alpha=100, lwd=3, add=F){
        #arguments passed directly to plotting functions

        #       
        if(!add){
                plot(self$time, exp(self$lq+log(self$N)),
                        type='l',
                        lwd=lwd,
                        col=col
                )
        } else{
                lines(self$time, exp(self$lq+log(self$N)),
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
plotRS = function(m=10^4, save=F){
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
        return(sam)
}
