tSource = function(fName){
        t=proc.time()
        source(fName)

        return( (proc.time()-t)/60 )
}
