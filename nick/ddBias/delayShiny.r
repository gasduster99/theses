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
        uniroot(function(FF){ eval(dFBdFexpr) }, c(0, 10))$root
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
        if( (t-aS)<1){
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
startWho = "./modsDDExpT45N300A0-1AS10K0.1/fit3KA_xi2.467_zeta0.538.rda"
dat = readRDS(startWho)
#aS = dat$aS
M = dat$M
dat$xi = NULL
dat$zeta = NULL

#
ui = dashboardPage(
	dashboardHeader(),
	dashboardSidebar(
		#Contrast
		sliderInput("sliderChi", "chi", min=0, max=1, step=0.05, value=1),
		#Growth & Maturity
		sliderInput("sliderAS","aS", min=0.1, max=10., step=0.1, value=dat$aS),
		sliderInput("sliderKappa","kappa", min=0.1, max=10, step=0.1, value=dat$kappa),
		#Recruitment
		sliderInput("sliderAlpha", "alpha", min=dat$M+0.1, max=exp(dat$lalpha)*1.5, step=0.05, value=exp(dat$lalpha)),
		sliderInput("sliderBeta" , "beta" , min=0, max=exp(dat$lbeta)*1.5, step=10^(floor(log10(exp(dat$lbeta)))-1), value=exp(dat$lbeta)),
		sliderInput("sliderGamma", "gamma", min=-2, max=2, step=0.10001, value=dat$gamma)
		
	),
	dashboardBody(
		fluidRow(column(12, plotOutput('bionum'))),
		fluidRow(column(12, plotOutput('nonlinear'))),
		fluidRow(column(12, plotOutput('catchsize'))),
		fluidRow(column(12,  plotOutput('refPlot')))
	)
)

#
server = function(input, output, session){
	#
	reactiveDat = reactive({
		#aS = input$sliderA
		
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
		dat$catch = fContrast(con)
                #
		dat$iterate()
		
		#return
		dat
	})	
	#
	output$bionum = renderPlot({
		#	
		dat = reactiveDat()
		layout(t(1:2))
		dat$plotQuan( function(B){B}, main="Biomass", ylim=c(0,max(dat$B)), xlab="Time", ylab="Biomass")
		dat$plotQuan( function(N){N}, main="Numbers", ylim=c(0,max(dat$N)), xlab="Time", ylab="Numbers")
	})
	#
	output$nonlinear = renderPlot({
		#
		dat = reactiveDat()
		#
		layout(t(1:2))
		#
		curve(SRR(x, dat), 0, 3*dat$B0, lwd=3, xlab="Biomass", ylab="Recruitment", main="Stock-Recruitment", n=1000)
		abline(0, dat$M, col='red')
		#
		ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0)
		curve(vbGrow(x, dat$kappa, dat$WW, dat$a0), 0, 15, lwd=3, xlab="Age", ylab="Biomass", main="VB Growth", ylim=c(0, dat$WW))
                segments(dat$aS, 0, dat$aS, ww)
		segments(0, ww, dat$aS, ww)
		points(dat$aS, ww, pch=19)
	})
	#
        output$catchsize = renderPlot({
                #       
                dat = reactiveDat()
                layout(cbind(c(1,2), 3))
		#
		#ww = vbGrow(dat$aS, dat$kappa, dat$WW, dat$a0) #WW*(1-exp(-kappa*a0))
		#dat$FMsy = FMsy(dat$M, dat$kappa, ww, dat$WW, exp(dat$lalpha), exp(dat$lbeta), dat$gamma)
		#
		par(mar=c(4,5,2,3))
		dat$plotQuan( function(catch, FMsy){catch*FMsy}, main="Fishing", ylim=c(0,1.5), xlab="Time", ylab="F")
		dat$plotQuan( function(B, catch, FMsy){B*catch*FMsy}, ylim=c(0,max(dat$B*dat$catch*dat$FMsy)*1.1), xlab="Time", ylab="Catch") 
        	par(mar=c(5,5,5,3))
		dat$plotQuan( function(B, N){B/N}, main="Average Size", ylim=c(0,max(dat$B/dat$N)), xlab="Time", ylab="Biomass Per Individual")
	})
	#
	output$refPlot = renderPlot({
		#
		dat = reactiveDat()
		##
		#RPs = NULL
		#RPs = rbind(RPs, c(dat$FMsy/dat$M, dat$zeta))
		#
		greys = rev(81-round(logseq(80, 1, 50))) #seq(60, 2, -1)
		howMany = length(dat$xi)	
		nWhite = max(howMany-length(greys)+1, 0)
		#
		plot(dat$xi, dat$zeta, pch=19, col=c(rep("white", nWhite), sprintf("grey%d",greys), "black"), xlab="Fmsy/M", ylab="Bmsy/B0", main="Reference Points")
		points(dat$xi[howMany], dat$zeta[howMany], pch=19, col='black')
		points(dat$xi[howMany], dat$zeta[howMany], col='red')
	})
}
#NOTE: Add varied contrast
#NOTE: plot R

#
#source('delayShiny.r'); shinyApp(ui, server)

#run app on host machine as below.
#connect via a web browser on LAN with host's IP:5050 in address bar
#source('delayShiny.r'); runApp(shinyApp(ui, server), host="0.0.0.0", port=5050)
