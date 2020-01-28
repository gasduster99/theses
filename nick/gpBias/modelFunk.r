#
dLikes = list(
        #
        LN = function(self, data){
                dnorm(log(data), log(self$q)+log(self$N), self$sdo, log=T)
        },
        #
        N = function(self, data){
                dnorm(data, self$q*self$N, self$sdo, log=T)
        }
)

#
qLikes = list(
        #
        LN = function(self, prob){
                qlnorm(prob, log(self$q)+log(self$N), self$sdo)
        },
        #
        N = function(self, prob){
                qnorm(prob, self$q*self$N, self$sdo)
        }
)
