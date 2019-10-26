organ = function(mixedList){
	out = list()
        for (j in seq(1, max(k))){out[[j]]=numeric(0)}
	for(ind in seq(1, length(mixedList))){
		jHat = length(mixedList[[ind]])
		out[[jHat]] = rbind(out[[jHat]], mixedList[[ind]])
	}
	return(out)
}
