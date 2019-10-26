## Active Learning for for regression by expected 
## improvement using Particle Learning on a simple 2-d
## exponential function

## load the plgp library
library(plgp)
library(tgp)
library(akima)

## close down old graphics windows
graphics.off()

## set up generation of Ys and Xs
runlock.basic = function(X)
{
  Y=numeric(nrow(X))
  for(i in 1:nrow(X)) {	
    write(c(6,X[i,]),file="input.tst",ncol=1)
    system("./RunLock input.tst output.tst 1")
    output = scan(file="output.tst")
##    Y[i] = output[1] + 10000*abs(output[2])+10000*abs(output[3])
    Y[i] = output[1] + 40000*(abs(output[2])>0) + 40000*(abs(output[3])>0) + 2*abs(output[2]) + 2*abs(output[3])
  }
  return(Y)
}
#runlock.basic = function(X) { 
#  if(is.matrix(X)) return(apply(X,1,sum)) 
#  else if(is.vector(X)) return(sum(X))
#  else { print("error in runlock"); return(-1) } }
#runlock.basic = function(X)
#{
#  Y=numeric(nrow(X))
#  for(i in 1:nrow(X)) {	
#    Y[i] = sum(X[i,])
#  }
#  return(Y)
#}
formals(data.GP.improv)$f <- runlock.basic

## set bounding rectangle for araptive sampling
rect <- matrix(c(rep(0,6),rep(20000,6)),ncol=2)
## rect <- matrix(c(rep(0,3),rep(20000,3)),ncol=2)
formals(data.GP.improv)$rect <- rect

## set a default prior
prior <- prior.GP(6)
prior$grate <- 50
## prior$drate <- -1/0.01
prior$bZero <- TRUE
formals(data.GP.improv)$prior <- prior
#formals(data.GP.improv)$oracle <- FALSE

## use a small LHS candidate set
formals(data.GP.improv)$cands <- 200

## set up start and end times
start <- 301
end <- 500

## do the particle learning
out <- PL(data=data.GP.improv, ## adaptive design PL via EI
          start=start, end=end, ## P=1,
          init=draw.GP,  ## init with Metropolis-Hastings
          lpredprob=lpredprob.GP, propagate=propagate.GP,
          prior=prior, addpall=addpall.GP, params=params.GP, cont=T)

## design a grid of predictive locations
XX <- dopt.gp(200, Xcand=lhs(200*10, rect))$XX
XXs <- rectscale(XX, rect)

## sample from the particle posterior predictive distribution
outp <- papply(XX=XXs, fun=pred.GP, quants=TRUE, prior=prior)

## extract the mean of the particles predictive
m <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) m <- m + outp[[i]]$m
m <- m / length(outp)

## extract the quantiles of the particles predictive
q2 <- q1 <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) {
  q1 <- q1 + outp[[i]]$q1
  q2 <- q2 + outp[[i]]$q2
}
q1 <- q1 / length(outp)
q2 <- q2 / length(outp)

## unscale the data locations
X <- rectunscale(pall$X, rect)

## plot the summary stats of the predictive distribution
par(mfrow=c(1,2))
image(interp.loess(XX[,1], XX[,2], m))
points(X)
image(interp.loess(XX[,1], XX[,2], q2-q1))
points(X)

## look at a historgram of the parameters
params <- params.GP()
Dev.new()
par(mfrow=c(1,2)) ## two plots
hist(params$d)
hist(params$g)

## plot EI progress
dev.new()
par(mfrow=c(1,2)) ## two plots
## plot the sampled points over time
plot(X[,1], main="sampled points",
     xlab="t", ylab="x1 & x2")
abline(v=start, col=3, lty=3)
lines((start+1):end, psave$xstar[,1])
points(X[,2], col=2, pch=18)
lines((start+1):end, psave$xstar[,2], col=2, lty=2)
legend("topright", c("x1", "x2"), col=1:2, pch=c(21,18))
## plot the max log ei over time
plot((start+1):end, psave$max.as, type="l", xlab="t",
     ylab="max log EI", main="progress meter")





## ## plot the summary stats of the predictive distribution
## pdf("nanprob_avg.pdf", width=5.5, height=5.5)
par(mfrow=c(2,3))
image(interp.loess(XX[,1], XX[,2], m),
      main="pred mean", xlab="x1", ylab="x2")
points(X[,1],X[,2],pch=19,cex=0.3)
image(interp.loess(XX[,3], XX[,4], m),
      main="pred mean", xlab="x3", ylab="x4")
points(X[,3],X[,4],pch=19,cex=0.3)
image(interp.loess(XX[,5], XX[,6], m),
      main="pred mean", xlab="x5", ylab="x6")
points(X[,5],X[,6],pch=19,cex=0.3)
image(interp.loess(XX[,1], XX[,4], m),
      main="pred mean", xlab="x1", ylab="x4")
points(X[,1],X[,4],pch=19,cex=0.3)
image(interp.loess(XX[,2], XX[,5], m),
      main="pred mean", xlab="x2", ylab="x5")
points(X[,2],X[,5],pch=19,cex=0.3)
image(interp.loess(XX[,3], XX[,6], m),
      main="pred mean", xlab="x3", ylab="x6")
points(X[,3],X[,6],pch=19,cex=0.3)
