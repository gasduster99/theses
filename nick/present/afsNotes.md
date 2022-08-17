---
title: Metamodeling for Bias Estimation of Biological Reference Points.
author: Nick Grunloh
documentclass: extarticle
geometry: margin=2cm
fontsize: 17pt
header-includes:
    - \usepackage{setspace}
    - \onehalfspacing
---

# Introduction

* Hello. My name is Nick Grunloh.
* Talking about: A Metamodeling approach for assessing RP bias in two parameter models of productivity.
<!--
* Talking about: A Metamodeling approach for assessing estimation bias in population dynamics models.
-->
* Collaboration with UC Santa Cruz, SWFSC, and funded by NOAA Sea Grant.

\clearpage

# Basic Modeling Structure

* Context: Single Species Surplus-Production Models.
* Production models are an admittadly simple setting, but...
	* have plenty of dark secrets that we dont tend to talk about. 
* Even being simple:
	* they capture many relevant dynamics for management sake
	* and are plenty instructive.

* General Structure:
	* Observe an index of abundance 
	* Assume the index is proportional to biomass with proportionality constant $q$.
	* $I_t$ forms a response variable with lognormal residuals.

* Most of the action here comes from the biomass process model.
	* Biomass is modeled as a (typically) nonlinear ODE.
	* Growth via a nonlinear production function, $P(B)$
	* Removals via natural mortality and catch.
		* Instantaneous removal rates lumped here under $Z(t)$. 

* For management mostly interested in Biological RP Inference.
	* I should say: RP are functions of productivity parameters.
* Commonly RPs are ways of noticing MSY.
	* Here I focus on two:
		* Fmsy: fishing rate to result in MSY (Relative. Fmsy/M)
		* Bmsy: population biomass at MSY (Relative. Bmsy/B0)

# RP Constraints 

* Conceptually $\frac{F^*}{M}$ and $\frac{B^*}{B_0}$ coexist in an entire 2D space.
* (Mangel et.al., 2013) Canadian Journal of Fisheries
	* Two parameter BH model: RP space is limited to a 1D curve
	* **Right:** Plot Relative Bmsy against Relative Fmsy
                
		* black: posterior samples of the RPs for a 3 parameter Shepherd-like model. (cowcod)
		* red: posterior samples of the RPs for a 2 param BH model.
		* the red posterior is squashed into the curve $\frac{1}{x +2}$
		* Refer to this subspace of RP's as the BH line.

	* **Next:** Mangel et. al. suggests looking into 3-parameter curves

# Schnute (1985)

* I'm working with a 3 parameter production function as developed in Schnute (1985)
	* A number of important 2 parameter special cases (Logistic, Ricker)
		* Most importantly here the Beverton-Holt when $\gamma=-1$.
* Generate a species off of the BH line, with the 3 parameter Schnute model.
* Simulate Schnute data and fit those data under the BH Model.
        * Observe how RPs under BH model are biased relative to true Schnute RPs.
* **Right Panel:** On the right you can see the Schnute production function and how it uses its third $\gamma$ parameter to get dramatically different productivity behavior.
* In terms of RPs these different behaviours move us off the BH line.

# \color{blue} Breadcrumb Slide

* Understanding the mapping of broad RP space onto these constrained 2 parameter spaces is complicated even in simple cases.
	* Chaos in the Dynamical System $~$ -Time Integrator Inaccuracy
	* Model Identifiability $~~~~~~~~~~~~~~~~$ -Global Optimization
* Production models are simplified places which are easier to hunt down the many computational issues, and are simple enough to make it possible to understand the mechanisms 
* At the link provided here you can see our anlysis of the mechanisms of Bias for the Schaefer Model.



\clearpage
# Simulation Design

* Again, we use Schnute to generate data broadly in RP-space.

* In order to do this, one would need to invert the relationship between RPs and productivity parameters.
* Schnute and Richards (1998) show that it's not analytically possible to invert this relationship
<!--* It's not analytically possible to invert the relationship between RPs and productivity parameters...--> 
and numerically inverting is unstable.
* However Schnute and Richards (1998) do provide some results which we have used to generate approximate LHS designs in RP space.
<!--partially write productivity parameters in terms of RPs.
Which can be used to generate approximate LHS designs broadly in RP space.
-->

* **Next:** With a design in place Schnute data can be generated for example in the upper right.
* Fitting the BH model against those data will neccessarily land somewhere on the $1/(x+2)$ BH line.  (Say the red dot)
* The aim is then to understand the behaviour of these bias arrow broadly in RP space.


* **Next:** A GP metamodel of the biases over RP space is fit using the observations at each design location
* Particular BH fits are only as helpful as their standard errors allow, but when you observe trends in RP bias on repeated sampling the metamodel can discern patterns of inferential bias and how it changes across RP space.

# Catch

* Assume synthetic catch series, To complete the model specification.
* I show Catch in red and the population in black.
* Two cases: Low and High contrast
* In all cases, the intitial biomass is fixed to $K$.
* Low Contrast:

	* Catch held at to come to equilibrium at MSY
	* low contrast, relatively low information setting
	* Exponential decay from K to Bmsy

* High Contrast:

	* fishing increases accelerates as technology and fishing techniques improve rapidly until management practices are applied to bring the stock into equilibrium at MSY.
	* high contrast, relatively high information setting
	* wiggles about until coming to equilibrium


<!--
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
-->


<!--
* Generate a species off of the Shaefer line, with a 3 parameter PT model.
* Same model structure as Shaefer, but production function is a 3 parameter curve.
* P even has a similar form.

        * PT has $\gamma$ parameter specifically to move the peak of the curve left or right.
        * moves the relative Bmsy off of the shaefer line at 1/2.

* Simulate PT data and fit those data under the Schaefer Model.

        * Observe how RPs under Schaefer model are biased relative to the true PT RPs.


* Extend simulation: replacing 3parm-PT with the 3parm-Dersio production function

        * Deriso has a number of 2 parameter special cases (including the Logistic)

                * Importantly Beverton-Holt and Ricker
-->
<!--
* Here I formulate a slightly simpler setting to investigate how a 2 parameter SRRs can bias RP estimation when fit to data generated under a 3 parameter production function.

* In particular I will generate data from a three parameter PT Production Model
        * PT is a 3 parameter generalization of the 2 parameter schaeffer model.
-->

<!--
* We can't observe all the fish in the sea, but thru the index
* We assume we observe biomass upto a proportionality constant $q$
* $q$ relates our index of abundance to actual biomass in the population.
* $I_t$ forms a response variable with lognormal residuals.

* We observe some index of population biomass up to a proportionality constant, $q$.
* Naturally the nonnegative index of abundance is observed with some uncertainty, which is typically assumed to have lognormal errors.

* Most of the action in these models comes in through a process model on Biomass.
* Biomass is modeled as a nonlinear ODE.
        * the population grows through a (typically non-linear) production function, P(B), and decreases as biomass is removed due to catch, C(t).
-->
<!--
Although this is an admittadly simple setting, production models have much 
to teach us and capture much of the realavant dynamics for the sake of management. 

* The modeling context here are single species population dynamics models as they might be used for managing fisheries.
* The simplest of these models (which really captures the essence of the managing objectives) is the surplus-production model.

* **Left Panel:** Data for a typical surplus-production model comes in the form of an index of abundance through time.
* **Right Panel:** The index is often observed alongside a variety of other known quantities, but at a minimum, each observed index will be observed in the presence of some known catch for the period.

* We can't observe all of the fish in the sea, but we can measure indicators of population biomass up to a proportionality constant, $q$.
* $q$ is the proportionality nuisance parameter (often called catchability) which relates our index of abundance to actual biomass in the population.
* And naturally the nonnegative index of abundance is observed with some uncertainty, which is typically assumed to have lognormal errors.

* Most of the action in these models comes in through a process model on Biomass.
* Biomass is modeled as a nonlinear ODE.
        * the population grows through a (typically non-linear) production function, P(B), and decreases as biomass is removed due to catch, C(t).
        * Production in this setting is defined as the net change in biomass due to basically all reproduction, maturation, and mortality processes other than the recorded fishing from humans.
        * Map the current biomass to some change (growth) in biomass

# Biological Reference points

* Reference points are simplified heuristic measures of population behavior, that are used to make decisions about how to manage the fishery.
* We want to mangage fisheries to allow (and promote) future productivity.
* The key idea is that we want to fish in a way to move the stable equilibrium of the population to a place along this curve that maximizes productivity in the steady state over time.
* Ex) Maximize simple yield at a particular moment V. Maximize sustainable yeild.

* The most common RPs are different ways of noticing that the population is at MSY.
* Any quantity decorated with a star reresent that quantity at MSY.
* Here I focus on the reference points Fmsy (fishing rate to result in  MSY) and Bmsy (biomass of the populaiton at MSY) (Or rather Bmsy as a fraction of K aka. Depletion at MSY)
-->

