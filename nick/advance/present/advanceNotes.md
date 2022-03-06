---
title: Metamodeling for Bias Estimation of Biological Reference Points.
author: Nick Grunloh
documentclass: extarticle
geometry: margin=2cm
fontsize: 14pt
---

# Introduction

* Hello. My name is Nick Grunloh. 
* Thanks to everyone in attendance and thanks to my committee for coming together on $\pi$ day 2022 to listen to my talk today.
* I'll be talking to you today about a Metamodeling approach for assessing estimation bias in important population dynamics models (the Schaeffer model), but also some extensions to this work as I finish up my PhD.
* This is work in collaboration with NOAA NMFS, and largely funded by NOAA Sea Grant, and I'd like to thank my collaborators with NOAA.

# Data and Basic Modeling Structure

* The modeling context here is that of population dynamics model as they might be used for stock assessment.
* I'll be focusing on the 


# (Mangel et.al., 2013) Canadian Journal of Fisheries and Aquatic Science

* I'll drop us into that work in the setting of a BH-SRR production model  
	
	* The primary insight of Mangel et.al. that I want to focus on is that while the space of $\frac{B^*}{B_0}$ and $\frac{F^*}{M}$ RPs is an entire 2D space, we quickly limit ourselves as we model.
	* With a 2 parameter production function, such as BH, if M is fixed this RP space is limited to a 1D curve.
	* Further if steepness ($h$) is specificed the space is further reduced to a 0D point along that curve, and we may have unintentionally selected our RPs before model even meets data.

\clearpage

* **Right Panel:** On the right I show an emperical way of noticing this 
<!--
, by fitting a model using a 2 parameter SRR and a model using a 3 parametes SRR on Cowcod Rockfish data.
-->	
	* The black are posterior samples of the RPs for a 3 parameter Shepherd-like model fit to cowcod rockfish data
	* The red are posterior samples of the RPs for a BH model with fixed M

* **Next:**
	
	* Notice how the posterior has been squashed into the red curve $\left(\frac{1}{\frac{F^*}{M} +2}\right)$
	* Mangel et. al. suggests looking into 3-parameter curves to avoid back ourselves into a corner and unintentionally overspecifying our models in the prior.

# Pella-Tomlinson Production Model

* Here I formulate a slightly simpler setting for investigating how models in the 2 parameter limited setting might be biased when fit to data from a 3 parameter production function.
* In particular I have a PT Production Model
* The model starts with a nuisance observation layer
* Biomass is driven by the given ODE.
* Production is given by R; That's the PT 3 parameter production function.

* **Right Panel:** On the right you can see the PT production function and how it uses its third $\gamma$ parameter to lean the production function to the left or right.

	* When $\gamma=2$; PT=Logistic Prodcution function and the curve is symmetric about $B^*$.
	
* **Next:** Recall the logistic production function is paramterized in terms of:

	* the slope at the origin ($r$) as seen in the slope of the blue line
	* and the right-hand x-intercept (K) as seen with the purple vertical line

<!--
* RP can easily be seen as a function of these parameters

	* $F^*=\frac{r}{2}$: slope of the red line
	* $B^*=\frac{K}{2}$: as seen in the vertical green line
-->

* **Next:** Due to the symmetry of the shaefer model, the RP space is limited to the horizontal line 1/2 @ all $F^*$
	
	* For brevity, later I'll refer to this line in this space as the shaefer line.
	* PT can get off of this line by changing the $\gamma$ parameter to lean left or right.

# Simulation

* The goal is to investigate bias induced my fitting PT data with the restricted shaefer model.
* I simulate data off of this shaefer line and subsequently fit those data with a limited 2 paramter model <!--production function; that limited model is the shaefer model.-->
* Here I show a grid of location in RP space where I will simulate data and once the shaefer model is fit, that estimate will neccessarily have to fit on that $\frac{1}{2}$ line.
* In particular I'll point out these 4 red X's in the 4 corners, Since I will refer back to these examples in about 2 slides.
* these are examples of large model misspecification.

# Catch

* I've looked at this basic setup across a range of different fishing behaviors.

	* Here I show 3 different fishing patterns
	* I've parameterized fishing relative to $F^*$ so that this first constant line at 1 on the left represents constant fishing at $F^*$.
	* The second curve represents a more typical ramp up of fishing and subsequent backing off toward $F^*$.
	* And the third curve is probably a completely uncommon fishing behavior just to see.

* **Next** For Brevity here we will only look at results of the constant fishing case.
* We will see that the detail you can get out of this simulation setting is plenty rich and the simplicty of this constant catch is helpful for understanding the mechanisms by which bias is induced when fitting a two parameter production function to even slightly more complicated data.

# Curves

* This slide visualizes the posterior fit of those data in the 4 large model misspecification corners of RP space
* In all cases the red lines represent the posterior fit of the Scheafer model, and the black is the truth of each of quantity.
* On the left I show the production curves of those 4 corner max misspecification fits with each plot in the relative position of each corner.
* In the middle column of small plots I shoe the posterior fit to biomass
* In the far right column of small plots I show the posterior fit to depletion

* **Left:** So starting on the left its immediatly possible to notice a few trends from these fits of the prodction function.
	
	* When the data are generated with $\frac{B^*}{B_0}>1/2$ (above the scheaffer line) $F^*$ is over estimated
		* We see that in the slope at the origin being too steep
	* When the data are generated with $\frac{B^*}{B_0}>1/2$ (above the scheaffer line) $F^*$ tends to be under estimated
	* Looking at these pictures you can start to understand why
		* When fishing is held at $F^*$ the population simply declines exponentially from $K$ to $B^*$.
		* The model only observes the right half of the true SRR
		* Due to the leaning of the true PT curves, and the symmetry of the logistic parabola, the logistic curve is learning about its slope at the origin entirely from data where depletion$>\frac{1}{2}$, and above the schaefer line PT is steeper that on the right half than it is on the left, and so we over estimate $F^*$ for data generated above the line.
		* The vice versa phenomena occurs below the schaffer line.
		* Data is only obsevered on the right half of the production function $\Rightarrow$ PT is shallower on the right than on the left $\Rightarrow$ and so the logistic parabola estimate tends to under estimate $F^*$.
	* Thats my $F^*$ story, but we can also observe some trends in biomass RPs.
	* Notice that the fits tend to match up the location of the humps fairly well
	* So the model tend to be estimating $B^*$ fairly well, but since we are fitting a restrictive parabola something has to give and you can see that we are totally missing on $K$.
* **Biomass:** In the center column we can see for the most part we doing ok on Biomass
* **Depletion:** But when you rescale things to consider depletion, we can see that the posterior estimate of ratio of biomass realtive to $K$, is ending up completely wrong for the highly misspecified cases.
* So that can be a small lesson that even if our models manage to predict faily well, model misspecification can completely lead us to incorrect inferences for some quantites. 

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

# Summary

* The statistician George Box famously said "All model are wrong, but some models are useful"
* My question is how useful are our models?, and if our models are useful. How useful are they? and when are they useful?
* This is a simulation based method for starting to understand that those questions.
* We want to expand these methods to more commonly used production functions
	* BH and Ricker
* and we want to try to get a similar analysis of these assumptiomtions when embeded inside of an age structured or delay difference model.
* but in this simple PT/Shaeffer setting we can see
	
	* as model misspecification increases biases in some quantities can become very large
	* estimates in $B^*$ are tending to be less sensative to model misspecification than $K$
	* and $F^*$ bias is going to tend to be strongly catch dependent


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




