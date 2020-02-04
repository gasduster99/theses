#
dLikes = list(
        #
        LN = function(self, data){
                dnorm(log(data), self$lq+log(self$N), self$sdo, log=T)
        },
        #
        N = function(self, data){
                dnorm(data, exp(self$lq)*self$N, self$sdo, log=T)
        }
)

#
qLikes = list(
        #
        LN = function(self, prob){
                qlnorm(prob, self$lq+log(self$N), self$sdo)
        },
        #
        N = function(self, prob){
                qnorm(prob, exp(self$lq)*self$N, self$sdo)
        }
)
