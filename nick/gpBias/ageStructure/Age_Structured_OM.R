###    AGE STRUCTURED OPERATING MODEL        ###
###            CHANTEL WETZEL                ###
###         chantell.wetzel@noaa.gov         ###


#Simulation Run Information=====================================================================

drive="C"
LH= "flatfish"
SRR = "BH" #"Shepard"
final.depl <- 0.25
tot.sims <- 1

#Set the directory==============================================================================

 #Create Files
 dir.create(paste(drive,":\\DBSRA\\",LH,sep=""))
 dir.create(paste(drive,":\\DBSRA\\",LH,"\\",model,sep=""))
 dir.create(paste(drive,":\\DBSRA\\",LH,"\\",model,"\\Run_Data",sep=""))

 wd <- paste(drive,":\\DBSRA\\",LH,"\\",model,sep="")
 setwd(wd)
 
 #Where to save output
 dynamics <- paste(wd,"\\Run_Data\\results_",depl.level,"_",sigma,sep="")
 #Load the biology 
 source(paste(drive,":\\DBSRA\\Code\\LH_parameter_values.R",sep=""))

#Read in the Catch History =======================================================================
 
 hist.catch<-c(rep(0,pre.fishery.yrs),rep(1000,fishery.yrs))  #dummy catch history 

#Paramaters========================================================================================
    
 NSIM <- as.vector(seq(1,tot.sims,1)) 
 
 #Bound for solving for R0
 first.R0 <- 10000
 low.bound <- 2000
 upper.bound <- 20000
 
  
 #add auto correlation in recruitment deviation
 rho=0.5 ; sigmaR <- 0.5
 
 #Survey Parameters
 Q<-1
 survey.time <- 0.5
 survey.CV <- 0.25
 
 #The pre.fishery years period is to spin the population out of equilibrium prior to fishing
 if (LH=="flatfish") {pre.fishery.yrs <- 30 }
 if (LH=="rockfish") {pre.fishery.yrs <- 100}
 fishery.yrs <- 50 
 setup.yrs <- pre.fishery.yrs + fishery.yrs
 total.yrs <- pre.fishery.yrs+fishery.yrs+1
 years <-c(seq(-pre.fishery.yrs,-1,1),seq(1,fishery.yrs+1,1))  
   
 
#Arrays----------------------------------------------------------------------------------------------------

  age.bins <- seq(1,num.ages,1)
  numbers <- array(0,dim=c(total.yrs,ages,sexes)) ; rownames(numbers)<-years
  biomass <- array(0,dim=c(total.yrs,ages,sexes)) 
  #Matrix Storage for output to R list
  like.matrix <- matrix(NA,length(NSIM),4)
  depl.hessian.run <- rep(NA,length(NSIM))
  SSB.matrix <- matrix(NA,total.yrs,length(NSIM))
  depl.matrix <- matrix(NA,total.yrs,length(NSIM))
  Ry.matrix <- matrix(NA,total.yrs,length(NSIM))
  rec.matrix <- matrix(NA,total.yrs,length(NSIM))
  vul.total.matrix <- matrix(NA,total.yrs,length(NSIM))
  ln.Ry.matrix <- matrix(NA,total.yrs,length(NSIM))
  index.expect.matrix <- matrix(NA,total.yrs,length(NSIM))
  vul.total.obs.matrix <- matrix(NA,total.yrs,length(NSIM))
  mid.phi.m.matrix <- matrix(0,length(len.step),ages)
  mid.phi.f.matrix <- matrix(0,length(len.step),ages)
  f.values.matrix <- matrix(0,(total.yrs-1),length(NSIM))
  fore.catch.matrix <- matrix(0, (total.yrs+1),length(NSIM))
  
#---------------------------------------------------------------------------------------------------
for (i in 1:length(NSIM))
 {

 nsim <- i
 
#Draw recruitment deviations--------------------------------------------------------------------------  
 set.seed(999+i)
 recdevs <- rnorm(total.yrs,0,sigmaR^2)
 autocorr <- rep(0,length(years))
  autocorr[1] <- recdevs[1]  
  for (e in 2:length(years)) { autocorr[e] <- rho*autocorr[e-1]+sqrt(1-rho*rho)*recdevs[e] }
 plot(years,autocorr,type='b',ylim=c(-0.35,0.35)); points(years,recdevs,col=2)
 
 set.seed(9999+i)
 survey.err <- rnorm(total.yrs,0,survey.CV^2)

#------------------------------------------------------------------------------------------------------
GetLen <- function()
{

 len<-matrix(NA,ages,sexes)
 mid.len<-matrix(0,ages,sexes)
 
 #L infinity (cm)
  Linf_f<-L1+((L2f-L1)/(1-exp(-kf*(a4-a3))))
  Linf_m<-L1+((L2m-L1)/(1-exp(-km*(a4-a3))))
  len.slope <- (L1-len.step[1])/a3
 
 #Length at the start of the year (cm)
  len[1:a.linear,]<-len.step[1]+len.slope*(seq(1,a.linear,1)-1)  #For smallest fish length is a linear function  
  #Growth based on the VB
  len[(a.linear+1):ages,2]<-Linf_m+(L1-Linf_m)*exp(-km*((seq(a.linear+1,ages,1)-1)-a3))
  len[(a.linear+1):ages,1]<-Linf_f+(L1-Linf_f)*exp(-kf*((seq(a.linear+1,ages,1)-1)-a3))
   
 #Mid-year lengths   (cm)
  mid.len[1:a.linear,]<-mean(len[1:2,1])+len.slope*(seq(1,a.linear,1)-1)   #For smallest fish length is a linear function 
  #Growth bases on VB
  mid.len[(a.linear+1):ages,2]<-len[(a.linear+1):ages,2]+(len[(a.linear+1):ages,2]-Linf_m)*(exp(-0.5*km)-1) 
  mid.len[(a.linear+1):ages,1]<-len[(a.linear+1):ages,1]+(len[(a.linear+1):ages,1]-Linf_f)*(exp(-0.5*kf)-1)
 
 output <- NULL
 output$len <- len
 output$mid.len<- mid.len
 return(output)
} 

#=================================================================================================================
Phi<-function() 
{ 
 mid.len <- GetLen()$mid.len
 len<-GetLen()$len 

 #St Dev of Mid-year lengths  (cm)
 mid.sigma<-matrix(0,ages,sexes) 
 sigma<-matrix(0,ages,sexes)
 phi.m<-matrix(0,length(len.step),ages) 
 phi.f<-matrix(0,length(len.step),ages) 
 mid.phi.m<-matrix(0,length(len.step),ages) 
 mid.phi.f<-matrix(0,length(len.step),ages)

 #CV about Length at the start of year
  sigma[1:a.linear,]<-len[1:a.linear,]*CV1
  sigma[(a.linear+1):floor(a4),] <-len[(a.linear+1):floor(a4),]*(CV1+(len[(a.linear+1):floor(a4),]-len[a.linear,])/(len[ceiling(a4),]-len[a.linear,])*(CV2-CV1))
  sigma[ceiling(a4):ages,]<-len[ceiling(a4):ages,]*CV2

 #CV about Length at the mid-year
  mid.sigma[1:a.linear,]<-mid.len[1:a.linear,]*CV1
  mid.sigma[(a.linear+1):floor(a4),] <-mid.len[(a.linear+1):floor(a4),]*(CV1+(mid.len[(a.linear+1):floor(a4),]-mid.len[a.linear,])/(mid.len[ceiling(a4),]-mid.len[a.linear,])*(CV2-CV1))
  mid.sigma[ceiling(a4):ages,]<-mid.len[ceiling(a4):ages,]*CV2

 #Start of year phi matrix
  for (b in 1:ages) 
   {   
    phi.m[1,b]<-pnorm(len.step[1+1],len[b,2],sigma[b,2])
    phi.f[1,b]<-pnorm(len.step[1+1],len[b,1],sigma[b,1])
    total.m <- phi.m[1,b]
    total.f<-phi.f[1,b] 
      for (a in 2:(length(len.step)-1)) 
      {
       p1 <- pnorm(len.step[a+1],len[b,2],sigma[b,2])
       p2 <- pnorm(len.step[a+1],len[b,1],sigma[b,1])
       phi.m[a,b] <- p1-total.m
       phi.f[a,b] <- p2-total.f
       total.m <- p1
       total.f <- p2
      } 
    phi.m[bin,b] <- 1-total.m
    phi.f[bin,b] <- 1-total.f
   }
        
  
 #Mid year phi matrix
  for (b in 1:ages) 
   {   
    mid.phi.m[1,b]<-pnorm(len.step[1+1],mid.len[b,2],mid.sigma[b,2])
    mid.phi.f[1,b]<-pnorm(len.step[1+1],mid.len[b,1],mid.sigma[b,1])
    total.m <- mid.phi.m[1,b]
    total.f<-mid.phi.f[1,b] 
     for (a in 2:(length(len.step)-1)) 
      {
       p1 <- pnorm(len.step[a+1],mid.len[b,2],mid.sigma[b,2])
       p2 <- pnorm(len.step[a+1],mid.len[b,1],mid.sigma[b,1])
       mid.phi.m[a,b] <- p1-total.m
       mid.phi.f[a,b] <- p2-total.f
       total.m <- p1
       total.f <- p2
      } 
    mid.phi.m[bin,b] <- 1-total.m
    mid.phi.f[bin,b] <- 1-total.f
   }
   
 output<-NULL
 output$mid.sigma<-mid.sigma
 output$sigma<-sigma
 output$phi.m<-phi.m
 output$phi.f<-phi.f 
 output$mid.phi.m<-mid.phi.m
 output$mid.phi.f<-mid.phi.f
 return(output)
}

#==========================================================================================
GetWght<-function()
{
 mid.phi.f<-Phi()$mid.phi.f
 mid.phi.m<-Phi()$mid.phi.m
 phi.f<-Phi()$phi.f
 phi.m<-Phi()$phi.m
  
 mid.wght<-matrix(0,ages,sexes)
 sample.mid.wght<-matrix(0,ages,sexes)
 wght<-matrix(NA,ages,sexes)
 wght.at.len<-matrix(0,length(len.step),sexes)
 mid.wght.at.len<-matrix(0,length(mid.len.step),sexes)

 #Virgin Weight @ Length (kg)
  wght.at.len[,2]<-(wght.coef.m*(len.step)^wght.exp.m)
  wght.at.len[,1]<-(wght.coef.f*(len.step)^wght.exp.f)
  mid.wght.at.len[,2]<-(wght.coef.m*(mid.len.step)^wght.exp.m)  #this is the value being outputted in ss3
  mid.wght.at.len[,1]<-(wght.coef.f*(mid.len.step)^wght.exp.f)

 #Virgin Weight @ Age  (kg)
  wght[,2]<-(t(phi.m))%*%mid.wght.at.len[,2]   
  wght[,1]<-(t(phi.f))%*%mid.wght.at.len[,1]
  mid.wght[,2]<-(t(mid.phi.m))%*%mid.wght.at.len[,2]
  mid.wght[,1]<-(t(mid.phi.f))%*%mid.wght.at.len[,1]    

 output <- NULL
 output$wght.at.len<-wght.at.len
 output$mid.wght.at.len<-mid.wght.at.len
 output$wght <- wght
 output$mid.wght <- mid.wght
 return(output)
} 

#================================================================================================
Fecundity<-function() 
{
 mid.wght <- GetWght()$mid.wght
 wght.at.len<-GetWght()$wght.at.len
 mid.wght.at.len<-GetWght()$mid.wght.at.len
 phi.f<-Phi()$phi.f 
 mid.phi.f<-Phi()$mid.phi.f

 mature.len<-1/(1+exp((ohm3)*(mid.len.step-ohm4)))#SS calcs using mid lengths
 mature.age<-(t(mid.phi.f))%*%mature.len  #SS uses mid phi values
 eggs<-ohm5+ohm6*wght.at.len[,1]  # eggs per kg

 fecund<-rep(0,ages)
 fecund<- (t(phi.f))%*%(mature.len*eggs*mid.wght.at.len[,1])

 output<-NULL
 output$mature.age<-mature.age
 output$mature.len<-mature.len
 output$fecund<-fecund 
 return(output)
}

#Virgin Population Structure========================================================================
Virgin <- function(R0) 
{ 
 fecund<-Fecundity()$fecund 
 len<-GetLen()$len 
 wght<-GetWght()$wght

 nnew<-matrix(NA,ages,sexes)
 nnew[1,]<-R0/2 
  
 for(a in 2:(ages-1))
  { 
   nnew[a,]<-nnew[a-1,]*exp(-m) 
   nnew[ages,]<-nnew[ages-1,]*exp(-m)/(1-exp(-m)) 
  } 
    
 #Virgin Biomass By Age
  bio<-nnew*wght  
  SSB0<-sum(nnew[,1]*fecund)
 
 output <- NULL
 output$biomass <- bio
 output$SSB0<-SSB0
 output$nnew<-nnew
 return(output) 
}

#Selectivity Function=============================================================
Selectivity<-function()
{
  mid.phi.m<-Phi()$mid.phi.m
  mid.phi.f<-Phi()$mid.phi.f
  selec <- matrix(NA,length(len.step),sexes)
  selec.age.m<-matrix(0,ages,1)
  selec.age.f<-matrix(0,ages,1) 
  
  #Double Normal Selectivity
  startbin <- 1
  peak <- fsp1
  upselex <- exp(fsp3)
  downselex <- exp(fsp4)
  final <- fsp6

  point1 <- 1/(1+exp(-fsp5)) 
  t1min <- exp(-((len.step[startbin]+1)-peak)^2/upselex)
  peak2 <- peak + 2 + (0.99*(len.step[length(len.step)]+1)-peak-2)/(1+exp(-fsp2))
  point2 <- 1/(1+exp(-final))
  t2min <- exp(-((len.step[length(len.step)]+1)-peak2)^2/downselex)
  t1 <- len.step+1-peak
  t2 <- len.step+1-peak2
  join1 <- 1/(1+exp(-(20/(1+abs(t1)))*t1))
  join2 <- 1/(1+exp(-(20/(1+abs(t2)))*t2))
  asc <- point1 +(1-point1)*(exp(-t1^2/upselex)-t1min)/(1-t1min)
  if (fsp5 <= -999) {asc =  exp(-(t1^2)/upselex)}
  dsc <- 1 +(point2-1)*(exp(-t2^2/downselex)-1)/(t2min-1)
  if (fsp6 <- -999) {dsc = exp(-(t2^2)/downselex)}

  selec[,1] <- asc*(1-join1)+join1*(1-join2+dsc*join2) 
  selec[,2]<-selec[,1]


 #Mid-year Selectivity by Age
  selec.age.m <-(t(mid.phi.m))%*%selec[,2]
  selec.age.f <-(t(mid.phi.f))%*%selec[,1]
  
output<-NULL
output$selec<-selec
output$selec.age.f<- selec.age.f
output$selec.age.m<- selec.age.m
return(output)
}  

#=========================================================================================================================================
Obs_Selectivity<-function()
{
 mid.phi.m<-Phi()$mid.phi.m
 mid.phi.f<-Phi()$mid.phi.f
 
 obs.selec<-matrix(0,length(mid.len.step),sexes) 
 obs.selec.age.m<-matrix(0,ages,1)
 obs.selec.age.f<-matrix(0,ages,1)

 #Survey Selectivity pattern Double Normal 
  startbin <- 1
  peak <- ssp1
  upselex <- exp(ssp3)
  downselex <- exp(ssp4)
  final <- ssp6

  point1 <- 1/(1+exp(-ssp5))
  t1min <- exp(-((len.step[startbin]+1)-peak)^2/upselex)
  peak2 <- peak + 2 + (0.99*(len.step[1]+1)-peak-2)^2/(1+exp(-ssp2))
  point2 <- 1/(1+exp(-final))
  t2min <- exp(-((len.step[length(len.step)]+1)-peak2)^2/downselex)
  t1 <- len.step+1-peak
  t2 <- len.step+1-peak2
  join1 <- 1/(1+exp(-(20/(1+abs(t1)))*t1))
  join2 <- 1/(1+exp(-(20/(1+abs(t2)))*t2))
  asc <- point1 +(1-point1)*(exp(-t1^2/upselex)-t1min)/(1-t1min)
  dsc <- 1 +(point2-1)*exp(-t2^2/downselex-1)/(t2min-1)

  obs.selec[,1] <- asc*(1-join1)+join1*(1-join2+dsc*join2) 
  obs.selec[,2]<-obs.selec[,1]
 
 #Mid-year Selectivity by Age
   obs.selec.age.m <-(t(mid.phi.m))%*%obs.selec[,2]
   obs.selec.age.f <-(t(mid.phi.f))%*%obs.selec[,1]

  output<-NULL
  output$obs.selec<-obs.selec
  output$obs.selec.age.m<-obs.selec.age.m
  output$obs.selec.age.f<-obs.selec.age.f
  return(output)
}  

#Recruits Spawning biomass  Vulnerable biomas  ===========================================================
Update_Dynamics <- function(R0, catch=hist.catch)
 {
 phi.f <-Phi()$phi.f
 phi.m <-Phi()$phi.m
 mid.wght<-GetWght()$mid.wght
 wght <- GetWght()$wght
 mid.wght.at.len<-GetWght()$mid.wght.at.len
 wght.at.len<-GetWght()$wght.at.len
 mid.len<-GetLen()$mid.len
 wght<-GetWght()$wght
 len<-GetLen()$len
 fecund<-Fecundity()$fecund
 mid.phi.m<-Phi()$mid.phi.m
 mid.phi.f<-Phi()$mid.phi.f
 sigma.len<-Phi()$sigma.len
 SSB0<-Virgin(R0)$SSB0
 selec<- Selectivity()$selec
 selec.age.f<- Selectivity()$selec.age.f
 selec.age.m<- Selectivity()$selec.age.m
 
 numbers[1,,] <- Virgin(R0)$nnew
 biomass[1,,] <- Virgin(R0)$biomass
 Ry<-matrix(0,total.yrs,1)
 Ry[1]<-R0/2   ;  rownames(Ry)<-years
 SSB<-matrix(0,total.yrs,1)
 SSB[1]<-SSB0  ;  rownames(SSB)<-years
 catch.at.age <- array(NA,dim=c(total.yrs,ages,sexes))  ; rownames(catch.at.age)<-years
 catch.at.len <- array(NA,dim=c(total.yrs,length(len.step),sexes))
 f.values <- rep(0,(total.yrs-1))
 catch.wght.values <- rep(0,(total.yrs-1)) 
 z.rate <- array(NA,dim=c(total.yrs,ages,sexes)) ; rownames(z.rate)<-years
 redo <- rep(0,total.yrs)

 for(y in 1:(pre.fishery.yrs+fishery.yrs))
 {  
      Findf <- function(f)
      {
        #Catch at age
        z.rate[y,,1] <- m+selec.age.f*f
        z.rate[y,,2] <- m+selec.age.m*f
        z.m <- (1-exp(-(z.rate[y,,2])))/(z.rate[y,,2])
        z.f <- (1-exp(-(z.rate[y,,1])))/(z.rate[y,,1])
        #Catch at Age
        catch.at.age[y,,1] <- max(0,f*(numbers[y,,1]*selec.age.f)*z.f)
        catch.at.age[y,,2] <- max(0,f*(numbers[y,,2]*selec.age.m)*z.m)
        #Catch At Length
        mid.temp.f <- numbers[y,,1]*z.f
        mid.temp.m <- numbers[y,,2]*z.m
        catch.at.len[y,,1] <- ((mid.phi.f*selec[,1])%*%(mid.temp.f))
        catch.at.len[y,,2] <- ((mid.phi.m*selec[,2])%*%(mid.temp.m))
        
        #Catch in Weight by Sex, mid.wght (41X2) calculated in the GetWght() function  
        catch.wght <- max(0,f*(sum(mid.wght.at.len[,1]*catch.at.len[y,,1])+sum(mid.wght.at.len[,2]*catch.at.len[y,,2])))    

        output <- NULL
        output$catch.at.age <- catch.at.age
        output$catch.at.len <- catch.at.len
        output$catch.wght <- catch.wght
        output$z.rate <- z.rate
        return(output)
      } #End FindF function 
      
      Obj.Fun.F <- function(f)
      {
        obj.fun.f <- (Findf(f)$catch.wght-hist.catch[y])^2
        return(obj.fun.f) 
      }
      f.value <- nlminb(0,Obj.Fun.F,lower=0,upper=4)
      f <- f.value$par 
  
      f.values[y] <- f  
      z.rate <- Findf(f)$z.rate
      catch.at.age <- Findf(f)$catch.at.age
      catch.at.len <- Findf(f)$catch.at.len
      catch.wght.values[y] <- Findf(f)$catch.wght
    
      #Update the numbers and remove the catch by applying the solved for f value
      for(a in 2:ages)
      {
       numbers[y+1,a,2] <- max(0,numbers[y,a-1,2]*exp(-m-selec.age.m[a-1]*f))
       numbers[y+1,a,1] <- max(0,numbers[y,a-1,1]*exp(-m-selec.age.f[a-1]*f))
      } 
       numbers[y+1,ages,2] <- max(0,numbers[y+1,ages,2] + numbers[y,ages,2]*exp(-m-selec.age.m[ages]*f)) 
       numbers[y+1,ages,1] <- max(0,numbers[y+1,ages,1] + numbers[y,ages,1]*exp(-m-selec.age.f[ages]*f)) 
    
      #Spawning biomass    
      SSB[y+1]<-0
      for (a in 2:ages)
      SSB[y+1]<-max(0,SSB[y+1]+numbers[y+1,a,1]*fecund[a])
      if (SSB[y+1] < 0) { SSB[y+1] <- 0 }

      if (SRR == "BH") {
      #Expected (and then realized) recruitment
      Ry[y+1]<-(4*h*(R0/2)*SSB[y+1])/(SSB0*(1-h)+SSB[y+1]*(5*h-1))
      Ry[y+1]<-Ry[y+1]*exp(-0.5*(sigmaR^2))*exp(recdevs[y+1]) }
      
      if (SRR == "Shepard") {
      #Expected (and then realized) recruitment
      Ry[y+1]<-#ADD SHEPARD'S SRR 
      Ry[y+1]<-Ry[y+1]*exp(-0.5*(sigmaR^2))*exp(recdevs[y+1]) }
         
      numbers[y+1,1,]<-Ry[y+1]
         
 } #closes yearly loop
 
 output <- NULL
 output$f.values <- f.values 
 output$z.rate <- z.rate
 output$catch.wght.values <- catch.wght.values 
 output$selec.age.m<-selec.age.m
 output$selec.age.f<-selec.age.f
 output$catch.at.age<-catch.at.age
 output$catch.at.len <- catch.at.len
 output$numbers <- numbers 
 output$SSB<-SSB
 output$Ry <- Ry
 output$redo <- redo
 return(output)
}

#=================================================================================================================================

 #Objective Function to solve for R0
 Test.R0 <- function(R0) {
  SSB<- Update_Dynamics(R0,catch=hist.catch)$SSB
  if ( SSB[pre.fishery.yrs+fishery.yrs+1]/SSB[pre.fishery.yrs] < final.depl)
    { obj.fun <- 10000000 }
    
  if ( SSB[pre.fishery.yrs+fishery.yrs+1]/SSB[pre.fishery.yrs] >= final.depl)
   {obj.fun <- (SSB[pre.fishery.yrs+fishery.yrs+1]-final.depl*SSB[pre.fishery.yrs])^2 }
  return(obj.fun)}

 #Find R0 that results in the correct depletion level
  r0.value=nlminb(first.R0,Test.R0,lower=low.bound,upper=upper.bound) 
  R0 <- r0.value$par

 # Create Virgin population
  dyn <- Update_Dynamics(R0,catch=hist.catch)
  virgin <- Virgin(R0)
  numbers[1,,] <- virgin$nnew
  biomass[1,,] <- virgin$biomass
  SSB0 <- virgin$SSB0 
  SSB<-dyn$SSB
  Ry <- dyn$Ry
  numbers <- dyn$numbers
  selec.age.f <- dyn$selec.age.f
  depl<-SSB/SSB[pre.fishery.yrs]


#===================================================================================================================================
DoSurvey <- function() 
{
 obs.selec.age.m<-Obs_Selectivity()$obs.selec.age.m
 obs.selec.age.f<-Obs_Selectivity()$obs.selec.age.f
 obs.selec<-Obs_Selectivity()$obs.selec
 catch.age.len <- dyn$catch.age.len
 numbers <- dyn$numbers
 catch.at.age<-dyn$catch.at.age
 selec <- Selectivity()$selec
 selec.age.f <- dyn$selec.age.f
 selec.age.m <- dyn$selec.age.m
 mid.phi.m<-Phi()$mid.phi.m
 mid.phi.f<-Phi()$mid.phi.f
 mid.wght.at.len<-GetWght()$mid.wght.at.len
 mid.wght <- GetWght()$mid.wght
 f.values <- dyn$f.values

 vul.bio.obs<-matrix(0,ages,1)
 vul.total.obs<-matrix(0,total.yrs,1)
 rownames(vul.total.obs)<-years
 index.expect<-matrix(0,total.yrs,1) ; rownames(index.expect)<-years
 survey.catch.age.len <-array(NA,c(total.yrs,ages,length(len.step),sexes))
 calc1 <- matrix(NA,length(len.step),ages)
 calc2 <- matrix(NA,length(len.step),ages)
 
 #Used in Catch @ Age and Length Calculation
 calc1 <- obs.selec[,1]*mid.phi.f
 calc2 <- obs.selec[,2]*mid.phi.m
 #Used in  Vulnerable Biomass Calculation
 temp2<-(t(mid.phi.m))%*%(mid.wght.at.len[,2]*obs.selec[,2])
 temp1<-(t(mid.phi.f))%*%(mid.wght.at.len[,1]*obs.selec[,1])

 for (y in 1:(pre.fishery.yrs+fishery.yrs))
  {
     vul.bio.obs <- (temp1)*numbers[y,,1]*exp(-survey.time*(m+selec.age.f*f.values[y])) +
                         (temp2)*numbers[y,,2]*exp(-survey.time*(m+selec.age.m*f.values[y]))      
     vul.total.obs[y] <- sum(vul.bio.obs[,])     
     #index.expect[y]<-Q*vul.total.obs[y]*exp(0.5*survey.CV-(survey.CV^2)/2)
     index.expect[y]<-Q*vul.total.obs[y]*exp(survey.err[y]-(survey.CV^2)/2)
     
    #Catch @ age and length
    for (a in 1:ages)
    {   
      survey.catch.age.len[y,a,,1] <- calc1[,a]*(numbers[y,a,1]*exp(-survey.time*(m+selec.age.f[a]*f.values[y])))      
      survey.catch.age.len[y,a,,2] <- calc2[,a]*(numbers[y,a,2]*exp(-survey.time*(m+selec.age.m[a]*f.values[y])))
    }        
  } 

 survey.data <- cbind(years,rep(1,total.yrs),rep(2,total.yrs),index.expect,rep(survey.CV,total.yrs))
 #survey.data <- survey.data[survey.duration,1:5]

 output <- NULL
 output$survey.catch.age.len <-survey.catch.age.len
 output$index.expect <- index.expect
 output$vul.bio.obs <- vul.bio.obs
 output$vul.total.obs <- vul.total.obs
 output$survey.data <- survey.data
 return(output)
}
#=====================================================================================================================

#Format Data: currently in the form for SS
 survey.data <- round(DoSurvey()$survey.data,2)

#=====================================================================================================================

#Store Values from OM population====================================================================================== 
#Dump the import data into storage matrices
  rec.matrix[,nsim] <- recdevs
  SSB.matrix[,nsim] <- dyn$SSB   
  depl.matrix[,nsim] <- dyn$depl  
  Ry.matrix[,nsim] <- 2*dyn$Ry
  ln.Ry.matrix[,nsim] <- log(2*dyn$Ry)
  vul.total.matrix[,nsim] <- dyn$vul.total
  index.expect.matrix[,nsim] <- DoSurvey()$index.expect
  vul.total.obs.matrix[,nsim] <- DoSurvey()$vul.total.obs
  mid.phi.m.matrix <- Phi()$mid.phi.m
  mid.phi.f.matrix <- Phi()$mid.phi.f
  f.values.matrix[,nsim] <- dyn$f.values
 
sim_info <- list()
sim_info[[1]] <- rec.matrix
sim_info[[2]] <- SSB.matrix
sim_info[[3]] <- depl.matrix
sim_info[[4]] <- Ry.matrix
sim_info[[5]] <- ln.Ry.matrix
sim_info[[6]] <- vul.total.matrix
sim_info[[7]] <- index.expect.matrix
sim_info[[8]] <- vul.total.obs.matrix
sim_info[[10]] <- mid.phi.f.matrix
sim_info[[11]] <- mid.phi.m.matrix
sim_info[[12]] <- Selectivity()$selec
sim_info[[13]] <- Obs_Selectivity()$obs.selec
sim_info[[14]] <- GetLen()$mid.len
sim_info[[15]] <- GetWght()$mid.wght
sim_info[[16]] <- GetWght()$mid.wght.at.len
sim_info[[17]] <- Fecundity()$fecund
sim_info[[18]] <- Fecundity()$mature.len
sim_info[[20]] <- dyn$selec.age.m
sim_info[[21]] <- dyn$selec.age.f
sim_info[[22]] <- obs.data$store.probs
sim_info[[23]] <- f.values.matrix
sim_info[[24]] <- like.matrix
sim_info[[25]] <- depl.hessian.run



save(sim_info,file=dynamics)

 
} #end simulation loop
