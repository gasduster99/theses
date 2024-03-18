rm(list=ls())

#https://stackoverflow.com/questions/53199645/r-an-interactive-graph-of-a-function-with-sliders-in-plotly
library(shiny)
library(shinydashboard)
#library(ggplot2)

#
library(plgp)
library(crch)
library(gMOIP)
library(Ryacas)
library(pracma)
library(stringr)
library(rootSolve)

#
source('ddClass0.1.1.r')

#
#DELAY FUNCTIONS
#

#https://www.wolframalpha.com/input?i=solve+w*a*B*%281-gamma*beta*B%29%5E%281%2Fgamma%29%2Bk*%28%28W*
#B = (1 - (((F + M) (F + k + M))/(a (F w + k W + M w)))^γ)/(β γ)
Bbar = ysym("(1 - (((F + M)*(F + k + M))/(alpha*(F*w + k*W + M*w)))^gamma)/(beta*gamma)")
Bbar_r = as_r( with_value(Bbar, "F", ysym("FF")) )
BBar = function(FF, M, k, w, W, alpha, beta, gamma, BbarX=Bbar_r){ eval(BbarX) }
#
FBbar = ysym("F*B")
FBbar = with_value(FBbar, "B", Bbar)
FBbar = with_value(FBbar, "F", ysym("FF"))
#
dFBdF = deriv(FBbar, "FF")
dFBdF_r = as_r(dFBdF)
FMsy = function(M, k, w, W, alpha, beta, gamma, dFBdFexpr=dFBdF_r){
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10), tol=.Machine$double.eps)$root
}
#
NBar = function(FF, M, k, w, W, alpha, beta, gamma){
        #
        BB = BBar(FF, M, kappa, ww, WW, alpha, beta, gamma)
        (alpha*BB*( 1-beta*gamma*BB )^(1/gamma))/(M+FF)
}
NCar = function(x, FF, M, k, w, W, alpha, beta, gamma){
        #
        #BB = BBar(FF, M, kappa, ww, WW, alpha, beta, gamma)
        (alpha*x*( 1-beta*gamma*x )^(1/gamma))/(M+FF)
}
#beta does not matter for either of these
#alpha|gamma, Fmsy
getAlphaFmsy = function(FF, M, k, w, W, beta, gamma, dFBdFexpr=dFBdF_r){
        #
        #capture.output(
        out <- tryCatch({
                uniroot(function(alpha){ eval(dFBdFexpr) }, c(eps(), 100))$root
        }, error=function(err){
                out = NA
        })#, file="/dev/null")
        #
        return(out)
}
#gamma|alpha, Fmsy
getGammaFmsy = function(FF, M, k, w, W, alpha, beta, dFBdFexpr=dFBdF_r){
        #FF = FMsy(M, k, w, W, alpha, beta, gamma)
        uniroot(function(gamma){ eval(dFBdFexpr) }, c(-10, 10))$root
}

#beta determines Bzero
getBeta = function(B0, M, k, w, W, alpha, gamma){
        f = function(b){ BBar(0, M, k, w, W, alpha, b, gamma) - B0 }
        uniroot(f, c(0, 10), tol=eps())$root
}

#gamma|alpha, zeta
getGammaZeta = function(zeta, FF, M, k, w, W, alpha, beta){
        f = function(g){
                BBar(FF, M, k, w, W, alpha, beta, g)/BBar(0, M, k, w, W, alpha, beta, g) - zeta
        }
        uniroot(f, c(-10, 10))$root
}

#
getZeta = function(FF, M, k, w, W, alpha, beta, gamma){
        BBar(FF, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
}

#
invert = function(zeta, xi, B0, M, k, w, W, alphaStart=1, betaStart=1, gammaStart=1){
        #
        alpha = alphaStart
        beta  = betaStart
        gamma = gammaStart
        #
        xTol = 10^-3
        zTol = xTol
        bTol = 10
        xiHat = xi+xTol*10
        zetaHat = zeta+zTol*10
        #
        while( abs(xiHat-xi)>=xTol | abs(zetaHat-zeta)>=zTol ){
                #
                gamma = getGammaZeta(zeta, xi*M, M, k, w, W, alpha, beta)
                alpha = getAlphaFmsy(xi*M, M, k, w, W, beta, gamma)
                beta  = getBeta(B0, M, k, w, W, alpha, gamma)
                #
                BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
                xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
                zetaHat = BBar(xiHat*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
                #
                #xis = c(xis, xiHat)
                #zetas = c(zetas, zetaHat)
                #i = i+1
        }
        #
        return( data.frame(alpha=alpha, beta=beta, gamma=gamma) )
}

#
vbGrow = function(a, k, W, a0){
        W*(1-exp(-k*(a-a0)))
}

#
getZetaBH = function(x, M, k, W, aS, a0){
        #
        w = vbGrow(aS, k, W, a0) #W*(1-exp(-k*a0))
        #
        gamma = -1
        alpha = getAlphaFmsy(x*M, M, k, w, W, 1, gamma)
        beta  = getBeta(B0, M, k, w, W, alpha, gamma)
        #
        BZero = BBar(0, M, k, w, W, alpha, beta, gamma)
        xiHat = FMsy(M, k, w, W, alpha, beta, gamma)/M
        zetaHat = BBar(x*M, M, k, w, W, alpha, beta, gamma)/BBar(0, M, k, w, W, alpha, beta, gamma)
        #
        return(zetaHat)
}
getZetaBH = Vectorize(getZetaBH, "x")

#DIFFERENTIAL EQUATION 

#
der = function(t, Y, lalpha, lbeta, gamma, aS, a0, WW, kappa, catch, B0){
        #linearly interpolate catches
        ft = floor(t)
        q  = (t-ft)
        Cl = catch[ft]
        Cu = catch[ft+1]
        C = q*Cu + (1-q)*Cl
        if(q==0){ C=Cl }
        #
        N = Y[1]
        B = Y[2]
        #
        if( (t-aS)<1){ #0.1
                Blag = B0
        }else{
                Blag = lagvalue(t-aS)[2]
        }
        #
        #R = exp(lalpha)*Blag*(1-exp(lbeta)*gamma*Blag)^(1/gamma)
        alpha = exp(lalpha)
        beta = exp(lbeta)
        R = alpha*Blag*(1-gamma*beta*Blag)^(1/gamma)
        #
        ww = vbGrow(aS, kappa, WW, a0) #WW*(1-exp(-kappa*a0))
        FF = C*FMsy(M, kappa, ww, WW, alpha, beta, gamma)
        #
        out = c(N=NA, B=NA)
        out[1] = R - (M+FF)*N
        out[2] = ww*R + kappa*(WW*N-B) - (M+FF)*B
        #
        return( list(out, ww*R, kappa*(WW*N-B)) )
}

#
SRR = function(B, mod){
        #
        alpha = exp(mod$lalpha)
        beta = exp(mod$lbeta)
        gamma = mod$gamma
        #
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
        #
        return(R)
}
SRR = Vectorize(SRR, 'B')

#
surplus = function(B, mod){
        #
        alpha = exp(mod$lalpha)
        beta = exp(mod$lbeta)
        gamma = mod$gamma
        #
	ww = vbGrow(mod$aS, mod$kappa, mod$WW, mod$a0)
	#
	FF = 0 #dat$FMsy
        R = alpha*B*(1-gamma*beta*B)^(1/gamma)
	N = R/(mod$M+FF)
	#Y = ww*R + mod$kappa*(mod$WW*N-B) - (M+FF)*B
	#Y = ww*R + mod$kappa*(mod$WW*R/(M+FF)-B) - (M+FF)*B
	Y = ww*R + mod$kappa*mod$WW*N - mod$kappa*B - mod$M*B - FF*B
        #fDat = function(x){ wwDat*SRR(x, dat) -dat$M*x-dat$kappa*x + 
	#	dat$kappa*dat$WW*NCar(x, 0, dat$M, dat$kappa, wwDat, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma) }
        # (alpha*x*( 1-beta*gamma*x )^(1/gamma))/(M+FF)
	return(Y)
}
surplus = Vectorize(surplus, 'B')

#CATCH

#
fContrast = function(con){
        TT = 31 #length(hake)
        tt = 1
        time = tt:TT
        #
        mid = round(TT/2)
        cMax = 1.6^con #2^con  
        cMin = 0.4^con #0.3^con #0.2^con
        b = log(cMin/cMax)/(tt-mid)
        a = exp(log(cMax)-b*mid)
        rSlope = (cMax-1)/(mid-1)
        rb = 2*rSlope*mid + cMax-rSlope*mid
        FtFmsy = (a*exp(b*time))*(time<=mid) + (-rSlope*time+rb)*(time>mid)
        #
        FtFmsy = c(FtFmsy, rep(1, 15))
}

fStep = function(){
	set.seed(1)
	FtFmsy = rep(0,46) #c(rep(5,15), rep(0, 16), rep(0, 15)) #+ rlnorm(46, 0, 0.25)
}

#
#SHINY
#

#https://stackoverflow.com/questions/53199645/r-an-interactive-graph-of-a-function-with-sliders-in-plotly
#https://mastering-shiny.org/basic-reactivity.html#simplifying-the-graph

#
#data <-  data.frame(x=c(1,2,3,4),y=c(10,11,12,13))
#startWho = "./modsDDExpT45N300A0-1AS10K0.1/fit3KA_xi3.493_zeta0.264.rda"
#startWho = "./modsDDExpT45N300A0-1AS10K0.1/fit3KA_xi3.225_zeta0.502.rda"
#startWho = "./modsDDExpT45N300A0-1AS10K0.1/fit3KA_xi2.478_zeta0.457.rda"
#path = "./modsDDExpT45N300A0-1AS10K0.1/"
#startWho = sprintf("%s/fit3KA_xi2.467_zeta0.538.rda", path)
path = "./modsDDExpT45N150A0-1AS2K0.1/"
startWho = sprintf("%s/fit3KA_xi1.463_zeta0.418.rda", path)
whoMight = as.data.frame(t(list.files(path=path, pattern=glob2rx("fit3KA*.rda"))))
colnames(whoMight)=whoMight
dat = readRDS(startWho)
dat$whoAmI = startWho
#aS = dat$aS
M = dat$M
dat$xi = NULL
dat$zeta = NULL

#
ui = dashboardPage(
	dashboardHeader(),
	dashboardSidebar(
		#Base Model NOTE: need to figure out how to update sliders to the new selected model
		#varSelectInput("who", "Base Model", whoMight, selected=startWho),
		
		#Contrast
		sliderInput("sliderChi", "chi", min=0, max=1, step=0.05, value=1),
		#Growth & Maturity
		sliderInput("sliderAS","aS", min=0.1, max=10., step=0.1, value=dat$aS),
		sliderInput("sliderKappa","kappa", min=0.1, max=5, step=0.1, value=dat$kappa),
		#Recruitment
		sliderInput("sliderAlpha", "alpha", min=dat$M+0.1, max=exp(dat$lalpha)*3, step=0.05, value=exp(dat$lalpha)),
		sliderInput("sliderBeta" , "beta" , min=0, max=exp(dat$lbeta)*1.5, step=10^(floor(log10(exp(dat$lbeta)))-1), value=exp(dat$lbeta)),
		sliderInput("sliderGamma", "gamma", min=-2, max=2, step=0.10001, value=dat$gamma)
		
	),
	dashboardBody(
		fluidRow(column(12, plotOutput('rowOne'))),
		fluidRow(column(12, plotOutput('rowTwo'))),
		fluidRow(column(12, plotOutput('rowThree'))),
		fluidRow(column(12, plotOutput('rowFour'))),
		fluidRow(column(12, plotOutput('rowFive')))
	)
)

#
server = function(input, output, session){
	#
	reactiveDat = reactive({
	
		##
		#if(input$who!=dat$whoAmI){
		#	whoNew = sprintf("%s/%s", path, as.character(input$who))
		#	dat = readRDS(whoNew)
		#	dat$whoAmI = whoNew 
		#}
			
                #Contrast
		con = input$sliderChi
		#Growth
		dat$aS = input$sliderAS
                dat$kappa = input$sliderKappa
                #Recruitment
		dat$alpha = input$sliderAlpha
		dat$lalpha= log(input$sliderAlpha)
		dat$beta  = input$sliderBeta
		dat$lbeta = log(input$sliderBeta)
		dat$gamma = input$sliderGamma	
		
		#
                ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0) #WW*(1-exp(-kappa*a0))
                dat$FMsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
		dat$xi   = rbind(dat$xi, dat$FMsy/dat$M)
		dat$zeta = rbind(dat$zeta, getZeta(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma))
		dat$catch = fStep() #fContrast(con)
                #
		dat$iterate()
		
		#return
		dat
	})	
	#
	output$rowOne = renderPlot({
		#	
		dat = reactiveDat()
		layout(t(1:2))
		dat$plotQuan( function(B){B}, main="Biomass", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass")
		dat$plotQuan( function(N){N}, main="Numbers", ylim=c(0,max(dat$N)), xlab="Time", ylab="Numbers")
	})
	#
	output$rowTwo = renderPlot({
		#
		dat = reactiveDat()	
		#
		layout(t(1:2))
		##
		#curve(SRR(x, dat), 0, 3*dat$B0, lwd=3, xlab="Biomass", ylab="Recruitment", main="Stock-Recruitment", n=1000)
		#abline(0, dat$M, col='red')
		#
		ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
		BMsy = BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)	
		g = function(x){BBar(x, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)}
		g = Vectorize(g, 'x')
		maxF = uniroot(g, c(dat$FMsy, exp(dat$lalpha)))$root
		#print(maxF)
		FFs = seq(0, maxF, length.out=1000)
		#f = function(x){surplus(x, dat)/surplus(BMsy, dat)*BMsy*dat$FMsy}
		#curve(f(x), 0, dat$B0, lwd=3, xlab="Biomass", ylab="Equilibrium Surplus Biomass", main="Yield Curve", n=1000)
		plot(g(FFs), FFs*g(FFs), type='l', lwd=3, xlab="Biomass", ylab="Equilibrium Surplus Biomass", main="Yield Curve")
		segments(BMsy, 0, BMsy, BMsy*dat$FMsy)
		points(BMsy, BMsy*dat$FMsy, pch=19)
		#abline(0, dat$FMsy)
		#curve(BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma))
		#rug( BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma), lwd=3 )
		#rug( BBar(0, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma), lwd=3 )
		#abline(h=BBar(dat$FMsy, dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)*dat$FMsy)
		#
		curve(vbGrow(x, dat$kappa, dat$WW, dat$a0), 0, 15, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, dat$WW))
                segments(dat$aS, 0, dat$aS, ww)
		segments(0, ww, dat$aS, ww)
		points(dat$aS, ww, pch=19)
	})
	#
        output$rowThree = renderPlot({
                #       
                dat = reactiveDat()
                #layout(cbind(c(1,2), 3))
		#
		#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0) #WW*(1-exp(-kappa*a0))
		#dat$FMsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
		#
		#par(mar=c(4,5,2,3))
		#dat$plotQuan( function(catch, FMsy){catch*FMsy}, main="Fishing", ylim=c(0,1.5), xlab="Time", ylab="F")
		#dat$plotQuan( function(B, catch, FMsy){B*catch*FMsy}, ylim=c(0,max(dat$B*dat$catch*dat$FMsy)*1.1), xlab="Time", ylab="Catch") 
        	#par(mar=c(5,5,5,3))
		#
		layout(t(1:2))
		#
		ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
		#
		curve(SRR(x, dat), 0, 3*dat$B0, lwd=3, xlab="Biomass", ylab="Recruitment #s", main="Stock-Recruitment", n=1000)
		abline(0, dat$M*(dat$M+dat$kappa)/dat$kappa/dat$WW/(1+dat$M*ww/dat$kappa/dat$WW), col='red')
		#abline(v=dat$B0)
		segments(dat$B0, 0, dat$B0, SRR(dat$B0, dat))
                #segments(0, SRR(dat$B0, dat), dat$B0, SRR(dat$B0, dat))
		points(dat$B0, SRR(dat$B0, dat), pch=19)
		#
		dat$plotQuan( function(B, N){B/N}, main="Average Size", ylim=c(0,max(dat$B/dat$N)), xlab="Time", ylab="Biomass Per Individual")
	})
	#
	output$rowFour = renderPlot({
		#
		dat = reactiveDat()
		#
		layout(t(1:2))
		dat$plotQuan( function(catch, FMsy){catch*FMsy}, main="Fishing", ylim=c(0,max(dat$catch*dat$FMsy)), xlab="Time", ylab="F")
                dat$plotQuan( function(B, catch, FMsy){B*catch*FMsy}, ylim=c(0,max(dat$B*dat$catch*dat$FMsy)*1.1), xlab="Time", ylab="Biomass", main="Catch")
	})
	#
	output$rowFive = renderPlot({
		#
		dat = reactiveDat()
		#
		howManyRP = length(dat$xi)
		howManyGrey = min(howManyRP, 50)
		greys = rev(81-round(logseq(80, 1, howManyGrey))) #seq(60, 2, -1)
		nWhite = max(howManyRP-howManyGrey+1, 0)
		#
		plot(dat$xi, dat$zeta, pch=19, col=c(rep("white", nWhite), sprintf("grey%d",greys), "black"), xlab="Fmsy/M", ylab="Bmsy/B0", main="Reference Points")
		points(dat$xi[howManyRP], dat$zeta[howManyRP], pch=19, col='black')
		points(dat$xi[howManyRP], dat$zeta[howManyRP], col='red')
	})
}
#NOTE: Add varied contrast
#NOTE: plot R

#
#source('delayShiny.r'); shinyApp(ui, server)

#run app on host machine as below.
#connect via a web browser on LAN with host's IP:5050 in address bar
#source('delayShiny.r'); runApp(shinyApp(ui, server), host="0.0.0.0", port=5050)


#dFBdF_stable = 
#expression((
#	1 - 
#	#(((FF + M) * (FF + k + M))/(alpha * (FF * w + k * W + M * w)))^gamma -
#	exp( gamma*( log(FF + M) + log(FF + k + M) )-(log(alpha)+log(FF * w + k * W + M * w)) ) -
#	(
#		#FF * (((FF + M) * (FF + k + M))/(alpha * (FF * w + k * W + M * w)))^(gamma - 1) * gamma * (alpha * (FF * w + k * W + M * w) * (FF + M + FF + k + M) - 
#		exp( log(FF) + (gamma - 1)*( log(FF+M)+log(FF+k+M)-log(alpha)-log(FF * w + k * W + M * w) ) + log(gamma) ) * 
#		#(alpha * (FF * w + k * W + M * w) * (FF + M + FF + k + M) -
#		exp( log(alpha) + log(FF * w + k * W + M * w) + log(2*FF+2*M+k) ) - 
#		#(FF + M) * (FF + k + M) * alpha * w)
#		exp( log(FF + M) + log(FF + k + M) + log(alpha) + log(w) )
#	) / (alpha * (FF * w + k * W + M * w))^2
#	) / (beta * gamma)
#)
#expression((
#	1 - 
#	(((FF + M) * (FF + k + M))/(alpha * (FF * w + k * W + M * w)))^gamma - 
#	(FF * (((FF + M) * (FF + k + M))/(alpha * (FF * w + k * W + M * w)))^(gamma - 1) * gamma * 
#	(alpha * (FF * w + k * W + M * w) * (FF + M + FF + k + M) - 
#	(FF + M) * (FF + k + M) * alpha * w)) / 
#	(alpha * (FF * w + k * W + M * w))^2
#	) / (beta * gamma)
#)
#NOTE: wrong
#expression(
#	1 - 
#	exp(gamma*( (log(FF+M)+log(FF+M+k)))-(log(alpha)+log(w)+log(FF+M+k*W/w))) - 
#	exp(
#		log(gamma)+log(FF)-log(alpha)-log(w) +
#		(gamma-1)*( (log(FF+M)+log(FF+M+k)))-(log(alpha)+log(w)+log(FF+M+k*W/w)) + 
#		log( 1 + (k*W/w)*(k-(k*W/w))/(FF+M+(k*W/w))^2 )
#	)
#)

