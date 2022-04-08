---
title: Metamodeling for Bias Estimation of Biological Reference Points.
author: Nick Grunloh
documentclass: extarticle
geometry: margin=2cm
fontsize: 17pt
header-includes:
    - \usepackage{setspace}
    - \singlespacing
---

# Introduction

* Hello. My name is Nick Grunloh. 
* Thanks for coming. $\pi$ day 2022!
* Collaboration with NOAA NMFS SWFSC, funded by NOAA Sea Grant.

# Outline 

* Metamodeling approach for assessing RP model misspecification in population dynamics models.
* Extensions to this work as I finish up my PhD.

<!--
* Thanks to everyone in attendance, and thanks to my committee for coming together today on $\pi$ day 2022 to listen to my advancement to candidacy talk.
* I'll be talking to you today about a Metamodeling approach for assessing estimation bias in population dynamics models, but also some extensions to this work as I finish up my PhD.
* This is work in collaboration with NOAA NMFS, and largely funded by NOAA Sea Grant, and I'd like to thank my collaborators with NOAA.
\clearpage
-->

\clearpage
# Data and Basic Modeling Structure

* Context:

	* single species population dynamics models
	* as might be used in managing fisheries

* Simplest model that captures the essence of the managing objectives is the Surplus-production model

* Data comes to us in two parts:

	* Two panels shown: example data from the classic Namibian Hake data set.
	* **Left** A measure of population biomass thru time (called an index of abundance) 
	* **Right** A Variety of other quantites, but at a minimum observe a matched time series of catches.

* We can't observe all the fish in the sea, but thru the index
* We assume we observe biomass upto a proportionality constant $q$
* $q$ relates our index of abundance to actual biomass in the population.
* $I_t$ forms a response variable with lognormal residuals.
 
<!--
* The modeling context here are single species population dynamics models as they might be used for managing fisheries.
* The simplest of these models (which really captures the essence of the managing objectives) is the surplus-production model.

* **Left Panel:** Data for a typical surplus-production model comes in the form of an index of abundance through time.
* **Right Panel:** The index is often observed alongside a variety of other known quantities, but at a minimum, each observed index will be observed in the presence of some known catch for the period.

* We can't observe all of the fish in the sea, but we can measure indicators of population biomass up to a proportionality constant, $q$.
* $q$ is the proportionality nuisance parameter (often called catchability) which relates our index of abundance to actual biomass in the population.
* And naturally the nonnegative index of abundance is observed with some uncertainty, which is typically assumed to have lognormal errors.
-->

* Most of the action in these models comes in through $B$ process.
* Biomass is modeled as a (typically) nonlinear ODE.
	* the population biomass grows via a production function, $P(B)$
	* biomass is removed due to catch, $C(t)$. 
	* Production maps the current biomass to some net growth of future biomass.
		* all Reproduction, maturation, and mortality processes other than fishing

# Schaefer Model

* Classic choice: $P(B)=$logistic growth curve => creates the Shaefer model.
* Parabola with parameters $(r, K)$
	
	* $r$ controls the maximum rate of productivity increase
	* $K$ is the so called "carrying capacity"

* Logistic models density dependence, and we can see the parameters on the parabola
	
	* low biomass: many resources & lack of competition => productivity increases quickly
	* biomass increases: more competition for resources => productivity decreases
	
* $K$ is a stable equilibrium point

	* population above $K$ decreases to return
	* population below $K$ increases to $K$

<!--
* One classic choice of the production function is the logistic growth curve.
* Choosing $P$ to be Logistic production function in the fisheries setting creates the Schaefer model. 
* In ecology logistic growth is very commonly parameterized in terms of the parameters $(r, K)$
* $r$ controls the maximum reproductive rate of the population in the absence of competition for resources (i.e. the slope at the origin). 
* $K$ is the so called "carrying capacity" of the population. 
* The quadratic shape of the logistic growth curve encodes density dependence
	* when the population biomass is low there is little competition for resources and growth is maximized (thus $r$).
	* when the population biomass is high there is competition for resources and growth declines until growth completely stops when the population reaches $K$.
	* In the absense of fishing $K$ is a stable equilibrium point, above $K$ production is negative to bring the population down, and below $K$ production is positive to bring the populaiton up to $K$.    
-->

# Biological Reference points

* We want to mangage fisheries to promote future productivity.
* $MSY$: Peak of production function.
* introduce fishing so as to move the stable equilibrium nearer to $MSY$

* **Next** 
\clearpage
* For management purposes we look at heuristic measures of the population called RP.
* Commonly RPs are ways of noticing MSY.
* Here I focus on two:

	* Fmsy: fishing rate to result in  MSY
	* Bmsy: biomass of the populaiton at MSY (argmax)
	* Relative Bmsy

* for a parabola, Relative Bmsy is just half of carrying capacity

<!--
* Reference points are simplified heuristic measures of population behavior, that are used to make decisions about how to manage the fishery.
* We want to mangage fisheries to allow (and promote) future productivity. 
* The key idea is that we want to fish in a way to move the stable equilibrium of the population to a place along this curve that maximizes productivity in the steady state over time. 
* Ex) Maximize simple yield at a particular moment V. Maximize sustainable yeild.
* The most common RPs are different ways of noticing that the population is at MSY.
* Any quantity decorated with a star reresent that quantity at MSY.
* Here I focus on the reference points Fmsy (fishing rate to result in  MSY) and Bmsy (biomass of the populaiton at MSY) (Or rather Bmsy as a fraction of K aka. Depletion at MSY)
-->
<!--
in time (and only for that moment) by fishing all available biomass in that moment. 
This strategy is penny-wise but pound-foolish (not to mention ecologically devastating) since it doesn’t leave biomass in the population to reproduce for future time periods. 
We seek to fish in a way that allows (or even encourages) future productivity in the population. 
This is accomplished by maximizing the equilibrium level of catch over time.
* Further if steepness ($h$) is specificed the space is further reduced to a 0D point along that curve, and we may have unintentionally selected our RPs before model even meets data.
, by fitting a model using a 2 parameter SRR and a model using a 3 parametes SRR on Cowcod Rockfish data.
* **Next:**
	
	* Notice how the posterior has been squashed into the red curve $\left(\frac{1}{\frac{F^*}{M} +2}\right)$
	* Mangel et. al. suggests looking into 3-parameter curves to avoid back ourselves into a corner and unintentionally overspecifying our models in the prior.
-->

# RP Constraints 

* Conceptually $F^*$ and $\frac{B^*}{B_0}$ coexist in an entire 2D space.
* (Mangel et.al., 2013) Canadian Journal of Fisheries and Aquatic Science

	* Two parameter BH model: RP space is limited to a 1D curve
	* **Right** Plot Relative Bmsy against Fmsy
		
		* black: posterior samples of the RPs for a 3 parameter Shepherd-like model. (cowcod)
		* red: posterior samples of the RPs for a 2 parameter BH model.
		* the red posterior is squashed into the curve $\frac{1}{x +2}$
	
	* Mangel et. al. suggests looking into 3-parameter curves


<!--
* The primary insight of Mangel et.al. that I'd like to focus on is that, if we are not careful, we quickly limit ourselves as we model.
* With a 2 parameter production function, such as BH, if M is fixed this RP space is limited to a 1D curve.
* **Right Panel:** On the right I show an emperical way of noticing this 
	
	* The black are posterior samples of the RPs for a 3 parameter Shepherd-like model fit to cowcod rockfish data
	* The red are posterior samples of the RPs for a BH model with fixed M
	* Notice how the posterior has been squashed into the red curve $\left(\frac{1}{\frac{Fmsy}{M} +2}\right)$
	* Mangel et. al. suggests looking into 3-parameter curves to avoid unintentionally overspecifying our models in the prior.
-->

* **Next:** 
\clearpage
* The Schaefer Model is a two parameter curve that suffers from a constrained RP-Space.

	* We already saw the constrain on the last slide
		
		* parobolic shape of logistic growth limits relative Bmsy to 1/2.
	
	* **Right**: This is the "Shaefer Line"
		
		* model misspecification in this context limits the space of RPs.
		* over constrains variance of estimates 
		* and induces severe bias in estimated RPs.
	
	* For managment define a specific (highly relavant) sense of model misspecification
		
		* Think of each point in this RP-space as a different species.
		* if spp not on line the model is misspecified.

# Outline

* In the next section I'll outline a simulaiton method for exploring the inferential effects of this type of model misspecification.

\clearpage
# Pella-Tomlinson Production Model

* Generate a species off of the Shaefer line, with a 3 parameter PT model.
* Same model structure as Shaefer, but production function is a 3 parameter curve.
* P even has a similar form.
	
	* PT has $\gamma$ parameter specifically to move the peak of the curve left or right.
	* moves the relative Bmsy off of the shaefer line at 1/2.

* Simulate PT data and fit those data under the Schaefer Model.

	* Observe how RPs under Schaefer model are biased relative to the true PT RPs.

<!--
* Here I formulate a setting for investigating how models in the 2 parameter limited setting might be biased when fit to data from a 3 parameter production function.
* In particular I have a PT Production Model
* The model starts with a nuisance observation layer
* Biomass is driven by the given ODE.
* Production is given by P; That's the PT 3 parameter production function.

* **Right Panel:** On the right you can see the PT production function and how it uses its third $\gamma$ parameter to lean the production function to the left or right.

	* When $\gamma=2$; PT=Logistic Prodcution function and the curve is symmetric about $Bmsy$.

* Simulate data under PT and fit those data under the Schaefer Model.
* Every point in RP-space corresponds to a set of parameters of the PT model.
-->

<!--
* RP can easily be seen as a function of these parameters

	* $F^*=\frac{r}{2}$: slope of the red line
	* $B^*=\frac{K}{2}$: as seen in the vertical green line


* **Next:** Due to the symmetry of the shaefer model, the RP space is limited to the horizontal line 1/2 @ all $F^*$
	
	* For brevity, later I'll refer to this line in this space as the shaefer line.
	* PT can get off of this line by changing the $\gamma$ parameter to lean left or right.

 at this basic setup across a range of different fishing behaviors.

	* Here I show 3 different fishing patterns
	* I've parameterized fishing relative to $F^*$ so that this first constant line at 1 on the left represents constant fishing at $F^*$.
	* The second curve represents a more typical ramp up of fishing and subsequent backing off toward $Fmsy$.
	* And the third curve is probably a completely uncommon fishing behavior just to see.

* **Next** For Brevity here we will only look at results of the constant fishing case.
* We will see that the detail you can get out of this simulation setting is plenty rich and the simplicty of this constant catch is helpful for understanding the mechanisms by which bias is induced when fitting a two parameter production function to even slightly more complicated data.
-->

\clearpage
# Catch

* Assume synthetic catch series, To complete the model specification.
* Data informs parameters differently than we normally expect.
* It is known: information content is about how biomass and catch series change.
	* not so much sample size
* I show Catch in red and the population in black.
* Two cases. All cases, Population initial condition starting at Carrying Capacity.
* Low Contrast:
	
	* Catch held at to come to equilibrium at MSY; Optimially managed species
	* low contrast, relatively low information setting
	* Exponential decay from K to Bmsy
	
* High Contrast:
	
	* fishing increases accelerates as technology and fishing techniques improve rapidly until management practices are applied to bring the stock into equilibrium at MSY.
	* high contrast, relatively high information setting
	* wiggles about until coming to equilibrium

	
<!--

* The information content of a given data series in this setting is known to be heavily dependent on catch and "contrast" in the time seeries.
* To control contrast I have catch parameterized in terms of a series of fishing rates relative to $Fmsy$.
* by varying $F/Fmsy$ you can get more or less contrast in the data.

* **Next:**

	* **Left:** By specifying $F/Fmsy=1$ it represents a low contrast, relatively low information, setting.
	* Below I show the catch in red and the biomass this induces in black.
	* We might call this a "one way trip"

	* **Right:** On the right $F/Fmsy$ is varied to induce a high contrast, relatively high information, setting.
	* This of the right panel as hypothetical stock where fishing rate accelerates as technology and fishing techniques improve rapidly until management practices are applied to bring the stock into equilibrium at MSY.
	* Below again you can see that this population is exposed to a varied catch history, which induce contrast in the generated indices and allows the fitting model to observe a decrease in the population followed by a rebuild of the stock
	* a so called "two-way trip" 

* I've looked other catch histories, but we'll focus on these tow extremes
-->
\clearpage
# Simulation

* Again, use PT to generate data broadly in RP-space.
* Shaefer model estimated RPs will neccessarily land somewhere on that $\frac{1}{2}$ Schaefer line. 
* Here I show my grid of RP locations for simulating data.
* **Next** 4 red X's in the corners are examples of large model misspecification of Schaefer
	
	* I show an example biomass series of in the high contrast setting in each corner.
	* $Bmsy/B0$: describes where the biomass comes to equilibrium
	* $Fmsy$: describes how quickly the stock responds to fishing and how fast it rebuilds.   

* **Next** At every simulation point, the parameters of the PT can be uniquely identified.

	* Analytical for PT
	* Generally inverting the RP-parameter relationship can be difficult.

\clearpage
# Metamodel

* Particular model fits are only as helpful as their standard errors allow, but when you observe trends in RP bias on repeated sampling broadly across RP space you can start to gain confidence in patterns of inferential bias.
* Here a Squared Exponential GP Metamodel is used as a flexible approximator for mapping the true PT RPs to the estimated RPs under limited setting of the Schaefer model. <!-- limitied estimate of RPs. parameter estimates when data are generated boardly over RP space.-->

* 

* Since the MLE estimator is a random variable it is important to propogate that uncertainty into the metamodel.
* The GP residual variation provides an ideal mechanism for propagating that uncertainty, to better understand biases.
* While the constrained RP space limits the extent of RP standard errors, accounting for estimate uncertainty has a smoothing effect, similar to a knugget, that improves the interpretability of RP biases. 

* 

* This metamodeling approach explicitly highlights the inferencial trade-offs imposed by productivity model misspecification in terms of specific bottom-line metrics that important to managing fisheries.
* While previous studies have considered the factors neccisary to estimate model parameters, the limiting constraints of RP model misspecification have not been explicitly considered, with a direct focus on RP estimation.
 

<!--
* The particular model fits from the last slide are only as helpful as their standard errors allow, but when you observe trends on repeated sampling you can start to gain confidence.
* With my grid of locations in RP space I am able to get a grid of fits, and across that grid I am able to establish bias tends in all of the major latent model quanities, and thats what I show on this slide.
-->

# Outline

* Split results into two parts
* First low contrast, low information, setting.

# Directionally

* Metamodel allows for a high level look at biases.
* What does inference do over the entire space of PT data
* Arrows indicate direction of mapping; color indicates the magnitude of bias.
* Give example. <!--, of data generated at location and mapped onto Shaefer line.-->

	* Map vertically to the horizontal Shaefer line
	* below line : underestimate $Fmsy$
	* above line: overestimate $Fmsy$

<!--
* Start with analysis of the low contrast, low information, catch setting. 

* Here I'm showing the entire space of PT data and when those data are fit with a Scheafer model how do RPs map onto the Schaeffer line?
* Give example, of data generated at location and mapped onto shaefer line.
* Arrows indicate direction of bias, and color indicates the magnitude of bias.

* Below the Shaefer line we underestimate $Fmsy$ and above the line we overestimate $Fmsy$. 
-->
\clearpage
# Components of Bias

* Split out thosse arrows into the components of bias.
* The color scale is in percent error.
* Red indicates over estimation, blue indicates underestimation by the Schaefer model.

* **Left** $Fmsy$ story
        
	* Below the line we underestimate $Fmsy$
	* Above the line we overestimate $Fmsy$
	* Larger Model misspecification => larger bias.

* **Right** Relative $Bmsy$ story

	* Law of Nature under Schaefer
	* Only looking at bias vertically
	* estimate must land on the line

<!--
* For all other plots Red indicates over estimation of the modeled quanity and blue indicates underestimates of the modeled quantity.
* **Left** In the bottom right you can see the generalized version of the our $F^*$ story
	* Below the line we underestimate $F^*$
	* Above the line we overestimate $F^*$
* **Right** I show the bias in $\frac{B^*}{B_0}$
	* this picture is a law of nature for this simulation setting (estimate must land on the line)
-->

# $Fmsy$ Curves

* How to understand $Fmsy$ bias?
* production functions in the 4 corners (examples of large model misspecification)

	* red lines represent the Schaefer fit
	* black represent the true PT production function
	* rug plot are the data

* Recall: low contrast, biomass is exponential decay from $K$ -> $Bmsy$.
\clearpage
* Tell top story

	* The model only observes the right portion of the true SRR	
	* Due to the leaning of the true PT curves, and the symmetry of the logistic 
		parabola, the logistic curve is learning about its slope at the origin 
		entirely from data where Biomass>Bmsy, and above the schaefer line the 
		PT is steeper that on the right half than it is on the left, and so we over 
		estimate $F^*$ for data generated above the line.  
	* The vice versa phenomena occurs below the schaffer line.  
	* $Fmsy$ is underestimated because of the shallow slopes on the right half of PT.

<!--
* This slide visualizes the posterior fit of those data in the 4 large model misspecification corners of RP space
* In all cases the red lines represent the posterior fit of the Scheafer model, and the black is the truth of each of quantity.
* I show the production curves of those 4 corner max misspecification fits with each plot in the relative position of each corner.

* When the data are generated with $\frac{B^*}{B_0}>1/2$ (above the scheaffer line) $F^*$ is over estimated
	* We see that in the slope at the origin being too steep
* When the data are generated with $\frac{B^*}{B_0}>1/2$ (above the scheaffer line) $F^*$ tends to be under estimated
* Looking at these pictures you can start to understand why
	* When fishing is held at $F^*$ the population simply declines exponentially from $K$ to $B^*$.
	* The model only observes the right half of the true SRR
	* Due to the leaning of the true PT curves, and the symmetry of the logistic parabola, the logistic curve is learning about its slope at the origin entirely from data where depletion$>\frac{1}{2}$, and above the schaefer line PT is steeper that on the right half than it is on the left, and so we over estimate $F^*$ for data generated above the line.
	* The vice versa phenomena occurs below the schaffer line.
	* Data is only obsevered on the right half of the production function $\Rightarrow$ PT is shallower on the right than on the left $\Rightarrow$ and so the logistic parabola estimate tends to under estimate $F^*$.
-->

# Ratio

* Recall: realtive $Bmsy$ pattern fixed in place
* Use Metamodel to split $Bmsy$ from carrying capacity (Constrained to divide back up).

	* $Bmsy$ shows large swaths of relatively little bias
	* Carrying capacity has substantial bias

		* Share one degree of freedom between the two qantities
		* K estimated to serve $Bmsy$
		* repeatedly observed pattern across all catch histories tried

<!--
* What's more interesting is what is the behavior of bias in the numerator ($B^*$) and the demoniator ($K$) not divided but independently
* Whatever the individual patterns are they need to divide back up to give this top middle picture.
* **Left:** On the left I show those quantities, the bias in the quantity $B^*$ is shown *top left*, and the bias in $K$ is shown on the *bottom left*.
* Interestingly, $B^*$ shows large swaths of relatively little bias, and most of the pattern in the top middle panel comes from bias in $K$.
* The world did not have to be this way, but this says that $Bmsy$ is often a robustly estimated quantity, but due to restrictive model misspecification, its a zero sum game and accuray in $B^*$ often come at the cost of estimates of $K$. 
-->
\clearpage
# Contrast

* I show $Fmsy$ and $Bmsy$
	
	* Constraint: carrying capacity bias is similar to previous slide

* When contrast is present, Bias is generally less over a large sets of largely misspecified settings.
	
	* As expected $Fmsy$ is largely estimated well in the presence of contrast.

**Next**

* The second row is the same high contrast setting, but appended an additionl period of catch near MSY.
	
	* Now, $Fmsy$ is more biased for a wide range of model misspecification
	* but, $Bmsy$ is better estimated now

* While contrast generally improves RP estimation, 

	* the over constrained Shaefer model lives in a zero sum world  
	 

<!--
# Heat Map

* The particular model fits from tyhe last slide are only as helpful as their standard errors allow, but when you observe trends on repeated sampling you can start to gain confidence.
* With my grid of locations in RP space I am able to get a grid of fits, and across that grid I am able to establish bias tends in all of the major latent model quanities, and thats what I show on this slide.
* **Top Right:** In the top right I'm showing the entire space of PT data and when those data are fit with a Scheafer model how do RP map onto the schaeffer line.
* For all other plots Red indicates over estimation of the modeled quanity and blue indicates underestimates of the modeled quantity.
* **Bottom Right** In the bottom right you can see the generalized version of the our $F^*$ story
	* Below the line we underestimate $F^*$
	* Above the line we overestimate $F^*$
* **Bottom Middle** Next to that we see a very similar pattern, but this time the quantity being graphed is $MSY$
\clearpage
* **Top Middle** I show the bias in $\frac{B^*}{B_0}$
	* this picture is a law of nature for this simulation setting (estimate must land on the line)
* Whats more interesting is whats is the behavior of bias in the numerator ($B^*$) and the demoniator ($K$) not divided but independently
* Whatever the individual patterns are they need to divide back up to give this top middle picture.
* **Left:** On the left I show those quantities, the bias in the quantity $B^*$ is shown *top left*, and the bias in $K$ is shown on the *bottom left*.
* Interestingly, $B^*$ shows large swaths of relatively little bias, and most of the pattern in the top middle panel comes from bias in $K$.
* The world did not have to be this way, but this says that $B^*$ is often a robustly estimated quantity, but due to restrictive model misspecification, its a zero sum game and accuray in $B^*$ often come at the cost of estimates of $K$. 
-->

\clearpage
# Summary

* "All model are wrong, but some models are useful"
* How useful is the Schaefer model?, and if useful, how useful? and when?
* This is a rich simulation based method understand that those questions at a deeper level.
	
	* The Shaefer model presents a simple setting to build the basic concepts at play, but 
	* I plan an extensions into a more broad class of production functions, as well as
	* extensions across models that include dynamics of individul growth and maturity.

* In this severly overconstrained settings we pay for our modeling mistakes primarily in estimate bias.
* In practice, when Schaefer model is unlikely to be correctly specified, one should at best expect to only reasonably estimate either $Bmsy$ or $Fmsy$ correctly depending on the particular degree of model misspecification.
* The observed contrast serves to distribute the available information among $B^*$ and $F^*$

	* So good models of catch are important to contextulize the interpretation of RP estimation.

<!--
* We want to expand these methods to more commonly used production functions
	* BH and Ricker
* and we want to try to get a similar analysis of these assumptiomtions when embeded inside of an age structured or delay difference model.
* but in this simple PT/Shaeffer setting we can see
	
	* as model misspecification increases biases in some quantities can become very large
	* estimates in $B^*$ are tending to be less sensative to model misspecification than $K$
	* and $F^*$ bias is going to tend to be strongly catch dependent
-->

# Outline: Start with the extension of growth and productivity

# Productivity Extension

* PT-Schaefer simulation setting explores a limited range of model misspecification matchups.

* Extend simulation: replacing 3parm-PT with the 3parm-Dersio production function

	* Deriso has a number of 2 parameter special cases (including the Logistic)
		
		* Importantly Beverton-Holt and Ricker

# Beverton-Holt

* Restricted inference under BH is of primary interest

	* Due to overwhelming popularity in Stock Assessment

* See 2 parmater form $\alpha$ & $\beta$.

	* $\alpha$ function similarly to $r$ in the Shaefer model
	* $\frac{1}{\beta}$ function similarly to the carrying capacity.

* **Right:** 

	* Again I show the constrained space=$\frac{1}{x+2}$
	* While equally constrained due to 2-parmater form, 
	* Cuts through RP-Space where we think many commonly assessed fish species exist in RP space. 

# Growth Extension

* Models that include individual growth and maturity can get more complex fast.
* Deriso-Schnute model: 
	* a relatively compact (fast) delay differential model to include growth effects, in the same simulation setting outlined here.
	* Models recruitment into numbers alongside a coupled biomass equation. 
	* Tracking numbers and biomass separately allows the model to account for the averge individual size and growth.

* Still RPs are largly determined by the recruitment function.

	* fills role of the production function 
	* very similar 2 v. 3 parameter parametric forms
	* Species Properties: $a0$, $\kappa$, $w_\infty$, $M$ estimated externally to the model.

* I will explore if the presence of these growth rate parameters effect RP bias

# Details

* Individuals simultaneously sexually mature, and start dieing and being fished at age $a_0$
* Individuals grow via von bertilanfy growth (Derived directly from).
* Biomass dynamics, can then be seen to describe biomass by an account 
	
	* B(new recruits), growth(existing biomass), mortality

# Outline 

* Given the importance of catch for interpreting RPs
* I have a catch model to better contextualize the dynamics model.

# Common Discretization

* The motivating biomass dynamics models are commonly discretized
* discretized for ease of implementation and data handling
* discrete version: catch observations assumed known and directly plugged in
* **Next**
        
	* The top requires instantaneous catch
	* The bottom only requires the aggregate catch over the year.
	* bottom matches the stucture of data as typiclly collected.
	* but the integral implies a useful latent stucture for modeling.
	
\clearpage
# Catch Interpolation

* **Right** I show the Namibian Hake catch again.
	
	* Each point represents the aggragate catch over a year.

* **Left** 
	* Model to get at the continuous latent structure of catch (consistent with statement of dynamics)
	* Continuous catch model is a linear spline.
	* Represents a sensible linear interpolation assumption on continuous catch.

* and so what of uncertainty in catch?
* Catch is commonly assumed known without uncertainty
* and there are biological hypotheses that could make use of jitters 
* A Statisticians eye will notice what looks like uncertainty in catch.
* so, I investigate what becomes of the uncertainty assumption?

# **Next**

* Working the linear spline assumption forward (integrating over time) to get a spline in terms of observables
<!-- * get another spline in terms of the observed quantities -->
* I show a model that allows uncertainty
* you can also turn uncertainty off (typically assumed)
* the plot below (**Left**) shows this model fit w/ and w/o assumpted uncertainty

	* Connect the dots vs. smooth

* **Next** (**Right**) 
	
	* Implied latent instantaneous catchs
	* catch uncertaitny smooths jitters and implies a sensible instantaneous catch.
	* assuming catch known implies wild oscillations on instantneous catch.

<!--
	* Ocillations neccessary on instantaneous scale to connect the dots on the aggregate scale (seen left) 
-->

* **Next** Back to the motivating models again, 
	
	* this exposes worrying assumptions
	* If we think that our discretizations bare any resemblance to the continuous motivating dynamics, then the assumption of catch known without uncertainty is also assuming that instantaneous catch can be negative (not an admissible assumption).  

* I'll investigate:

	* Models of various smoothness on instantaneous catch  	
	* how accounting for catch uncertainty affects inference in the biomass dynamics model 

# Timeline

<!-- * and further explore how assuming --> 
<!--
* Biomass dynamics models are commonly discretized, for ease of implementation and data handling
* The top is the motivating continuous formulation of dynamics
* The bottom is a common discretization
* The discrete version is susceptible to all sorts of numerical problems, but it can be convient
* **Next**

	* The top requires instaneous catch
	* The bottom only requires the aggregate catch over the year.
-->

<!--
* We are capable of estimating variance for catch in each year so, I've shown an option for how to handle uncertainty, and 
-->
<!--	
	* 
	* Catch is commonly assumed known without uncertainty
	* there are biological hypotheses
	* A Statisticians eye will notice what looks like uncertainty in catch.
-->

<!--
* **Figures**

	* **Left** shows the model on the observation level
	* **Right** shows the implied latent instantaneous catchs 

* The blue is the model fit with-out assumed uncertainty
* The black is the model fit with assumed uncertainty
* On the observation scale the blue simply connects the dots, but the implied instantaneous catch oscillated wildly
* Assuming catch uncertaitny smooths the observed jitters and implies a much more sensible instantaneous catch.
-->

<!--
* In the middle column of small plots I shoe the posterior fit to biomass
* In the far right column of small plots I show the posterior fit to depletion
	* Thats my $F^*$ story, but we can also observe some trends in biomass RPs.
	* Notice that the fits tend to match up the location of the humps fairly well
	* So the model tend to be estimating $B^*$ fairly well, but since we are fitting a restrictive parabola something has to give and you can see that we are totally missing on $K$.
* **Biomass:** In the center column we can see for the most part we doing ok on Biomass
* **Depletion:** But when you rescale things to consider depletion, we can see that the posterior estimate of ratio of biomass realtive to $K$, is ending up completely wrong for the highly misspecified cases.
* So that can be a small lesson that even if our models manage to predict faily well, model misspecification can completely lead us to incorrect inferences for some quantites. 


-->

<!--
# 78-82 Bars

* **Top Panel:** 
	
	* For each market category accounting for 99% of landings  
		* (blue) Proportion landings by weight
		* (red) Proportion samples by #

* **Bottom Panel:** Aggregated Species Compositions
	* Colors represent 13 select species (others grey)
	* Number above is the # species present
	* Hatching is MCAT nominal species 

* MCATs not pure
	* often nominal species is not even the major species
		* BCAC
		* BRWN

* Sampling Opportunistic
	* Sampling co-occurs with landings
	* Often as more species are present there are more samples
	* This is lucky for modeling
		* More samples than parameters (largely driven by spp)
		* Most landings are modeled (78-82: 96.8%)

* No sampling south of Conception

# 83-90 Bars

* Same picture but 83-90
	
	* (blue) Proportion landings by weight
	* (red) Proportion samples by #
	* Aggregated Species Compositions

* Top 99% of landings in more market categories
	
	* MCATs still largely impure

* 83-90: 98.3% of landings modeled


\clearpage

# Likelihood Forms

* First modeling choice: Pick a Likelihood

* Shelton etal. 2012 Fit Multinomial via the Multinomial-Poisson trans.
	
	* Piece together independent Poissons

* We are not limited to Multinomial distribution
	
	* quantify uncertainty (residual variability)
	* consider modeling overdispersion
	* additional parameter ($\phi$) to disentangle mean from variance

* $y_{ij}$: $i^{\text{th}}$ sample of the $j^{\text{th}}$ species' integer
weight

* Remove all other modeling decisions by modeling a single stratum
	
	* MCAT 250
	* Montery
	* Trawl
	* 1982/Q2

# Likelihood Graphs

* Fit models and look are how they predict

* **Left Pannel:** 95% HDI from each model along side observed sppComp data

	* black horizontal lines are observed species comps
	* blue is Possion (i.e. Multinomial) Model
	* red is Binomial
	* green is the Negative Binomial Model
	* yellow is the Beta-binomial Model

* **Right Panel:** Entire Beta-binomial predictive distribution

* Overdispersion is present (spp comps from [0,1])

* ~50 obsevations => 2.5 missing in 95% interval

	* Maybe NB missing a few to many, and BB missing a few to few
	* BB certaintly finding the most variance
	* split intervals but... very appropriate density

# Likelihood Table

* Consider MSE, DIC, WAIC, and Marginal Likelihood Bayesian Model Prob.

* Varied model selection criterion (Nothing is perfect!) 

* Consistent and large support for the Overdispersion Models
	
	* Most support for BB

* Moving forward I develop the BB model   

# Beta-Binomial Model

* A Full Operationalized Model!

* $y_{ijklm\eta}$: $i^{\text{th}}$ sample of the $j^{\text{th}}$ species' integer weight, in the $k^{\text{th}}$ port, caught with the $l^{\text{th}}$ gear, in the $\eta^{\text{th}}$ \mbox{quarter,} of year $m$, for a particular market \mbox{category.}

* Stratum $\mu$ linked to $\theta$ and observed cluster size ($n$) 

* Stratum $\sigma^2$ is largely a function of $\mu$ but with overdispersion $\rho$
	
	* $\rho\rightarrow0$: Binomial variance
	* $\rho\rightarrow1$: $n$ times Binomial variance

* Modeling of $\theta$ (all predictors are categorical):
	
	* Intercept
	* Additive offsets for: Species, Port, Gear
	* Consider multiple time models

\clearpage

# Time Models

* Bayesian Modeling 
	* Heirarchical v. Random Effect Disclaimer

* (M1) Fixed main effect time model
	
	* No pooling

* (M2) Random main effect time model
	
	* years/quarter pool separately

* (M3) Random main effects + random interaction

* (M4) Random interactions jointly pooled

* (M5) Random interactions quarterly variances pooling across years

* (M6) Random interactions yearly variances pooling across quarters 

# Priors

* Very diffuse priors

* Main effects diffuse Normals

* $\rho$ transformed to be a real number

	* $\text{logit}(\rho) \rightarrow (-3.91, 3.91)$
	* $\rho \rightarrow (0.02, 0.98)$  

* Any heirarchical variance gets the same IG prior
	
	* Considered others:
		* $\sqrt{v}~$~ Half-Cauchy$(10^{-2})$
		* $\sqrt{v}~$~ Unif$(0, 10^5)$

\clearpage

# Beta-Binomial Fits

* Fit model separately in 78-82 and 83-90 and compare model selection criterion

* 78-82:

	* Consistent support for more pooling
	* All measures point to (M4)

* 83-90:
	
	* Consistent support for interaction models
	* Uncertainty between (M3), (M4), and (M5)
	* Lesser support for (M6)

* We fit model (M4) everywhere
	
	* Stable and relatively fast model to fit
	* Given its support in 78-82, I am drawn to (M4)
		* Each time period seems to have a mind of its own  
\begin{eqnarray*}
&\beta^{(t)}_{m\eta} = \beta^{(y)}_{m} + \beta^{(q)}_{\eta} + \beta^{(y:q)}_{m\eta} & \\
&\beta^{(y)}_{m} \sim N(0, 32) & \\
&\beta^{(q)}_{\eta} \sim N(0, 32) & \\
&\beta^{(y:q)}_{m\eta} \sim N(0, v) &
\end{eqnarray*}

# ?? LUNCH ??

\clearpage

# Posterior Predictive Species Comps.

* Having settled on (M4) in both time periods, how do we build species comps?

* Inference results in samples from posterior distribution $P\Big(\mu_{jklm\eta}, \sigma^2_{jklm\eta} | y\Big)$

* Run samples back through BB likelihood to compute Monte Carlo integral and get posterior predictive distribution of sampled weight.

* Use draws from model posterior predictive weight to compute species comp. distribution
	
	* Plot shows average species compositions 
	* Full distribution for $y^*$ as well as $\pi^*$
	* Each sample sums to 1 and $\sum_j\mathbb{E}[\pi^*_j]$=1

* By adding an unobserved latent time period we can make out-of-sample predictions
	* (M4): unobserved $\beta^{(y)^*}$ and $\beta^{(q)^*}$

# Single Quarter Hindcast

* Recall for 1978-1982 there was no sampling south of point conception.

* Adding an unobserved year and quarter

	* make predictions for each species in each combo of: 
		* three observed gear groups 
		* three southern port complexes

\clearpage

# 78-82 Prediction

* Modeled MCATs

* MCATs in the rows **(ordered by landings)** w/ 3 nominal HDI prediction levels

	* For each stratum of each MCAT compare data to prediction intervals
	* Observed level should match Nomial
	* Prediction higher than nominal => Overfitting
	* Prediction lower => Underfitting (not enough residual variance) 
	
* Most do well
	* Average performance is reasonable
	* Note this is a unweighted, simple, average
	* More accurate would weight average by samples at each stratum 

* Particularly well in heavily landed stratum

	* correlation of sampling effort w/ landings

* Widow is a wild child
	* only example that is off by more than 5% points at any level

# 83-90 Prediction

* Same Table
	* Modeled MCATs **(ordered by landings)** w/ 3 nominal pred. levels

* Again most do well

* Recall landings were spread across more MCATs in 83-90
	* Enough samples to also model more MCATs

* Blackgill, Yellowtail, Cowcod: off by 5% points at some level
	* Negligable Landings

# Speciating Landings

* $\lambda_{\cdot klm\eta}$ is reported on landing reciepts

* $\lambda^*_{jklm\eta}$ stored in DB

* Aggregate to any level
	
	* across quarter, port complex, gear group
	* Also MCAT (I ran out of index variables :/)

* E.J. will show the speciated time series with predictive intervals 
	
	* summed across MCAT
	* as it might be used in asessment

\clearpage

# BMA Story

* Mentioned partial pooling thru time via heirarchical modeling

* But present system also pools in space
	* Given sparcity, it's entirely possible that we also need spatial pooling 

* I show MSE to demonstrate the biase/variance trade off
	
	* Pooling directly exchanges sample size (postior variance) with bias
	* **Far Right** Least Bias
	* **Far Left** Most Bias, but most data for small posterior variance
	* A practicle solution is somewhere in between

* [Bell number] Idea: Try all partions of port complexes

* $\text{B}_{10}=115975$
	
	* Too many
	* add Spatial Modeling Constraints 
	* Biogeography viewed through the lens of human behavior
		* sampling behavior
		* fisherman behavior

* $\bar{\text{B}}_{10}=61136$
	
	* Require partitions to be "small"
	* No super-grouping greater than 3 port complexes
	* points close in space behave similarly (smoothness)

* $\hat{\text{B}}_{10}=512$
	
	* Require continuous partitions (no leap frogging)
	* Like a GP continunity constraint

* $\hat{\bar{\text{B}}}_{10}=274$

	* Together we have "small" and "continuous" partitions
	* smoothness and continuity 
	* a computationally manageable set of models to compute

\clearpage

# BMA Math
	
* Defines a candidate model set

* We could just pick the single "best" model

	* defining "best" is hard
	* model selection criterion are imperfect 

* Prediction results are averaged results

# 78-82 BMA Results

* Describe plot

* Recall no sampling in the south
	* All latent structure filled in by predictive distribution in the south

* 250:
	* Marginal model probability 
		* N1: 32+14+13+12=71%
		* N2: 2+2+2+2=8%

* 253:
	* Central Block
	* BRG/BDG
	* Lump/Split CRS and ERK (among top 5 models; 58% model weight)
		* Split: 0.3448276
		* Lump: 0.6551724

* 269:
	* Sold on the BRG-BDG break

\clearpage

# 83-90 BMA Results

* Recall BRG-OSF Missing data

* 250:
	
	* Missing data story
		* Lump or
		* Quarantine	

* 956:
	
	* Lump/Split CRS and ERK
	* A break at Cape Mendicino
	* BRG/BDG/OSF Quarantine os missing data

* 269:
	
	* Piont Conception Break
	* Cape Mendicino Break


# Conclusions

* Using Bayesian models we have:
	
	* Account for overdispersion
	* Estimate uncertainty (full distribution)
	* Formal Mechanisms for pooling
	* provide structure for making out-of-sample prediction

* Future Modeling
	
	* Explore additional predictore in $\theta$
		* Landing weighting
		* Vessel Effects
		* Speceies:Gear interactions
	
	* Overdispersion Multivate models
		* Dirichelette-Multinomial Model
	
	* Maybe Time Series Models
	
	* Cluster and integrate out spatial parameters via DP?
-->

<!--
# 91-99 Bars

# 00-15 Bars

# $\rho$ Posterior

# $v$ Posterior

# Spp Comps Sum to 1
-->




