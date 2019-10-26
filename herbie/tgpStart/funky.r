rosenbrock = function(x){
        out = 100*(x[,1]^2 - x[,2])^2 + (x[,1] - 1)^2
        return(out)
}

rosenbrockHD = function(x){
        d = dim(x)[2]
        out = 0 
        for (i in seq(1, d-1)){
                inc = 100*(x[,i]^2 - x[,i+1])^2 + (x[,i] - 1)^2
                out = out + inc 
        }   
        return(out)
}

rastrigin = function(x){
        out = matrix(NaN, nrow=dim(x)[1], ncol=1)
        for (i in seq(1, dim(x)[1])){
                ex = x[i,]
                out[i] = 10*length(ex)+sum(ex^2-10*cos(2*pi*ex))
        }   
        return(out)
}

