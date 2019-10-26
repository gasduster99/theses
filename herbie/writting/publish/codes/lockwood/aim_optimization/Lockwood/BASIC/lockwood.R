runlock.basic = function(rates=c(rep(10000,6)))
{
  write(c(6,rates),file="input.tst",ncol=1)
  system("./RunLock input.tst output.tst 1")
  output = scan(file="output.tst")
  return(data.frame(cost=output[1],plume.a=output[2],plume.b=output[3]))
}




runlock.basic(c(353,6003,13526,2472,777,600))

runlock.basic(c(rep(100,6)))

