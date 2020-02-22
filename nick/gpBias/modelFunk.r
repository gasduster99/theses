#
dLikes = list(
        #
        LN = function(self, data){
                dnorm(log(data), self$lq+log(self$N), exp(self$lsdo), log=T)
        },
        #
        N = function(self, data){
                dnorm(data, exp(self$lq)*self$N, exp(self$lsdo), log=T)
        }
)

#
qLikes = list(
        #
        LN = function(self, prob){
                qlnorm(prob, self$lq+log(self$N), exp(self$lsdo))
        },
        #
        N = function(self, prob){
                qnorm(prob, exp(self$lq)*self$N, exp(self$lsdo))
        }
)
