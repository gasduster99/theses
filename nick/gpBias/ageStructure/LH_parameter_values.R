if (LH == "flatfish" ) {
 sexes=2
 ages <- 31 
 num.ages <- ages-1
 fishing.time <- -1 ; max.F <- 4 ; F.type <- 2 
 start.bias <- 1980         ;  full.bias <- 1990
 last.bias <- 2006          ;  last.no.bias <- 2009
 max.bias.adj <- 0.50
 pre.model.devs <- 1930

 len.step<-seq(from=10,to=64,by=2)      ; mid.len.step<-seq(from=11,to=65,by=2)  
 bin<-length(len.step)
         
 m <- 0.20                  ; h <- 0.875 
 wght.coef.m <- 0.00000717  ; wght.coef.f <- 0.00000342
 wght.exp.m <- 3.134        ; wght.exp.f <- 3.346
 L1 <- 24.6219              ; ohm3 <- -0.734  #slope of mat fnc
 L2f <- 55.4099             ; ohm4 <- 33.10   #length at 50% mat
 L2m <- 40.6664             ; ohm5 <- 1
 kf <- 0.143779             ; ohm6 <- 0
 km <- 0.299548             ; sigmaR <- 0.5         ; rho <- 0.50 #autocorrelation
 a3 <- 2.833                ; a4 <- 17.833               
 a.linear <- floor(a3)      ; len.slope <- 5.1613
 CV1 <- 0.08                ;  CV2 <- 0.08
 Amat <- 5
 
 Fmsy <- 0.2362918
 Fmsy.M = Fmsy/m
 Bpeak <- 0.25 
 
 fsp1 <- 43               ; ssp1 <- 33
 fsp2 <- 3                ; ssp2 <- 3 
 fsp3 <- 5.05             ; ssp3 <- 5.05  
 fsp4 <- 6                ; ssp4 <- 6
 fsp5 <- -12              ; ssp5 <- -12 
 fsp6 <- 70               ; ssp6 <- 70
 
 target <- 0.25
 
 #Forecast File Values=================================================================================================
 spr.target <- 0.25         ; bio.target <- 0.25
 
 selec.slope <- 33         ; selec.inflec <- 8
 survselec.slope <- 24      ; survselec.inflec <- 8
 
 ctl.rule.tgt<- 0.25        ; ctl.rule.thres <- 0.05
} 



#Rockfish Section========================================================================================
if (LH == "rockfish" ) { 
 sexes=2
 ages <- 101 
 num.ages <- ages-1
 fishing.time <- -1 ; max.F <- 4 ; F.type <- 2 
 start.bias <- 1980         ;  full.bias <- 1990
 last.bias <- 2006          ;  last.no.bias <- 2009
 max.bias.adj <- 0.50
 pre.model.devs <- 1930
 
 len.step<-seq(from=4,to=66,by=2)      ; mid.len.step<-seq(from=5,to=67,by=2)  
 bin<-length(len.step)
       

 m <- 0.05                  ; h <- 0.50 
 wght.coef.m <- 0.00000155  ; wght.coef.f <-0.00000155
 wght.exp.m <- 3.03         ; wght.exp.f <- 3.03
 L1 <- 6.64                 ; ohm3 <- -0.25  #slope of mat fnc
 L2f <- 59.844              ; ohm4 <- 40.5   #length at 50% mat
 L2m <-L2f-0.1346           ; ohm5 <- 1
 kf <- 0.1314               ; ohm6 <- 0
 km <- 0.257985             ; sigmaR <- 0.50    ; rho <- 0.50 #autocorrelation
 a3 <- 1                    ; a4 <- 80                  
 CV1 <- 0.08                ; CV2 <- 0.08
 Amat <- 9
 
 Fmsy <- 0.03186427
 Fmsy.M = Fmsy/m
 Bpeak <- 0.40            
 
 
 fsp1 <- 43               ; ssp1 <- 33
 fsp2 <- 3                ; ssp2 <- 3 
 fsp3 <- 5.05             ; ssp3 <- 5.05  
 fsp4 <- 6                ; ssp4 <- 6
 fsp5 <- -12              ; ssp5 <- -12 
 fsp6 <- 70               ; ssp6 <- 70
 
 target <- 0.40
 
 #Forecast File Values=================================================================================================
 spr.target <- 0.50         ; bio.target <- 0.40
 
 selec.slope <- 34          ; selec.inflec <- 4
 survselec.slope <- 34      ; survselec.inflec <- 4
 
 ctl.rule.tgt<- 0.40        ; ctl.rule.thres <- 0.10
 
}
