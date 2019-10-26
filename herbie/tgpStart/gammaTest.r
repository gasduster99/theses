rm(list=ls())

as = seq(0.001, 10, 0.5)
bs = seq(0.01, 10, 0.01)

flag=1
for (a in as){
	for(b in bs){
		print(a)
		print(b)
		print('')
		if(flag){
			curve(dgamma(x, a, scale=b), 0, 20, 1000,)
			flag=0
		}
		curve(dgamma(x, a, scale=b), 0, 20, 1000, add=T)
	}
}
