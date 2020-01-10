rm(list=ls())

load('save.RData')

#
wk = 0.001
gamma = -1
kappa = 0

#
h  = 0
sh = samples$sigma*(1-h)
W0 = ((1-sh)*wk + (1-kappa)*sh*wk)/(1-kappa*samples$sigma)
rh = (wk/W0)*(1-sh)/(1-h)
R0 = (1-(rh/samples$alpha)^gamma)*rh/(samples$beta*gamma) #(1-exp(gamma*(log(rh)-log(alpha))))*(rh/(beta*
B0 = (W0/wk)*(R0/(1-sh))
#
rStar = sqrt(samples$alpha*(1-samples$sigma))
hStar = rStar - (1-samples$sigma)
hStar = hStar/(1+hStar)
CStar = ((sqrt(samples$alpha) - sqrt(1-samples$sigma))^2) / (samples$beta)
RStar = ((sqrt(samples$alpha) - sqrt(1-samples$sigma))    / (samples$beta)) * sqrt(samples$alpha)
SStar =  (sqrt(samples$alpha) - sqrt(1-samples$sigma))    / (samples$beta*sqrt(1-samples$sigma))
BStar = SStar * (1+sqrt(samples$alpha*(1-samples$sigma))-(1-samples$sigma))

#
#Marc EJ Pic
#

#
Fmsy = -log(1-hStar)
EM = -log(samples$sigma) 
#
dev.new()
plot(Fmsy/EM, BStar/B0)

#
#R(t) v.S(t)
#

#
gamma = -1
S = 1:(max(postB))
R = matrix(NA, nrow=M, ncol=length(S))
Rtb = matrix(NA, nrow=length(S), ncol=2)
#
L = matrix(NA, nrow=M, ncol=length(S))
Ltb = matrix(NA, nrow=length(S), ncol=2)
#
D = matrix(NA, nrow=M, ncol=length(S))
Dtb = matrix(NA, nrow=length(S), ncol=2)
for(s in S){
	#
	R[,s] = samples$alpha*(s)*(1-samples$beta*gamma*s)^(1/gamma)
	Rtb[s,] = quantile(R[,s], c(0.025, 0.975))
	#
	#L[,s] = sqrt(samples$alpha*(1-samples$sigma))*s
	#L[,s] = rh*s
	#Ltb[s,] = quantile(L[,s], c(0.025, 0.975))
	##
	#D[,s] = R[,s]-L[,s]
	#Dtb[s,] = quantile(D[,s], c(0.025, 0.975))
}
#MSY = matrix(NA, nrow=M, ncol=1)
#Smsy = matrix(NA, nrow=M, ncol=1)
#Rmsy = matrix(NA, nrow=M, ncol=1)
#Bmsy = matrix(NA, nrow=M, ncol=1)
#for(m in 1:M){
#	MSY[m] = max(D[m,])
#	Smsy[m] = S[D[m,]==MSY[m]]
#	Rmsy[m] = R[D[m,]==MSY[m]]
#	Bmsy[m] = B[t] = Rmsy[m] + (((1-kappa)*wk + kappa*W[t-1])/W[t-1]) * samples$sigma*(Smsy)
#}

#
dev.new();
plot(S, colMeans(R), ylim=c(0, max(Rtb)), type='l')
polygon(c(S, rev(S)), c(Rtb[,1], rev(Rtb[,2])), border=NA, col=rgb(0,0,0,0.5))
#polygon(c(S, rev(S)), c(Ltb[,1], rev(Ltb[,2])), border=NA, col=rgb(0,0,0,0.5))
lines(S, colMeans(R), lwd=5)
#lines(S, colMeans(L), lwd=5)
##
#dev.new();
#plot(S, colMeans(D), ylim=c(0, max(MSY)), type='l')
#polygon(c(S, rev(S)), c(Dtb[,1], rev(Dtb[,2])), border=NA, col=rgb(0,0,0,0.5))
#lines(S, colMeans(D), lwd=5)
#points(Smsy, MSY)
##
#dev.new();
#hist(MSY)



