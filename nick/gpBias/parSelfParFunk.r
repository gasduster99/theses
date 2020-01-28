#
selfToPar = function(parNames){
        #check if variable names exist

        #
        parValues = c()
        for(pn in parNames){
                eval(parse( text=sprintf("parValues['%s']=self$%s", pn, pn) ))
        }
        #
        return(parValues)
}

#NOTE: parValues should be passed with names
parToSelf = function(parValues){
        #check is names exist

        #
        parNames = names(parValues)
        #
        for(pn in parNames){
                #update self
                eval(parse( text=sprintf("self$%s=parValues[pn]", pn) ))
        }
}
