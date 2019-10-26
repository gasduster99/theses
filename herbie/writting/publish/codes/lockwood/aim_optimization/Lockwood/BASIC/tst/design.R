## version to run on the Lockwood problem

## illustrates sequential design and optimization
## by active learning heuristics; the default is
## ALM, but you can also try ALC and EI

## load the libraries
library(dynaTree)
## dyn.load("../dynaTree/src/dynaTree.so")
library(tgp)

rect <- matrix(c(rep(0,6),rep(20000,6)),ncol=2)

## Initial (x,y) data
start <- 50
X <- dopt.gp(start, Xcand=lhs(10*start, rect))$XX
y <- numeric(start)
for(i in 1:start) {
  out <- runlock.basic(X[i,])
  y[i] <- out$cost + 10000*abs(out$plume.a)+10000*abs(out$plume.b)
}

## size of predictive grid and type of AS
ngrid <- 1000
method <- "ei" ## also try "alm", "alc", "ei", or "ieci"
prec <- 0.1 ## for ei
ieci <- alc <- ei <- FALSE
if(method == "alc") { alc <- TRUE
 } else if(method == "ei") { ei <- TRUE 
 } else if(method == "ieci") ieci <- TRUE
  
## PL fit to initial data
obj <- dynaTree(X=X, y=y)##, model="linear")

## determining the number of adaptive sampling rounds
end <- 1000

##
## Do the sequential design
##

track <- NULL

for(t in start:end){

  ## random predictive grid
  XX <- lhs(ngrid, rect)

  ## predict at the XX locations
  obj <- predict(obj, XX, quants=FALSE, ei=ei)
  if(alc) { obj <- alc(obj, XX, rect=rect)
  } else if(ieci) obj <- ieci(obj, XX)

  ## extract via ALM, ALC, EI-prec
  al <- alcalc(obj, method, prec)
  m <- which.max(al)
  track <- c(track, al[m])
  xstar <- drop(obj$XX[m,])
  out <- runlock.basic(xstar)
  ystar <- out$cost + 10000*abs(out$plume.a)+10000*abs(out$plume.b)

  ## plot the progress
#  plotprogress(obj, xstar, ystar, method, track, f1d)
  
  ## update the fit for the next round
  obj <- update(obj, matrix(xstar, nrow=1), ystar, verb=1)

  print(c(ystar, min(obj$y)))
}

## free the particle cloud on the C-side
## deletecloud(obj); obj$num <- NULL

##dev.off()
