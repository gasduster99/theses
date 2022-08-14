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
* Talking about: A Metamodeling approach for assessing estimation bias in population dynamics models.
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
* Commonly RPs are ways of noticing MSY.
	* Here I focus on two:
		* Fmsy: fishing rate to result in MSY (Relative Fmsy)
		* Bmsy: biomass of the populaiton at MSY (Relative Bmsy)

# RP Constraints 

* Conceptually $\frac{F^*}{M}$ and $\frac{B^*}{B_0}$ coexist in an entire 2D space.
* (Mangel et.al., 2013) Canadian Journal of Fisheries
	* Two parameter BH model: RP space is limited to a 1D curve
	* **Right:** Plot Relative Bmsy against Relative Fmsy
                
		* black: posterior samples of the RPs for a 3 parameter Shepherd-like model. (cowcod)
		* red: posterior samples of the RPs for a 2 param BH model.
		* the red posterior is squashed into the curve $\frac{1}{x +2}$

	* **Next:** Mangel et. al. suggests looking into 3-parameter curves

# \color{blue} Breadcrumb Slide

* Understanding the mapping of broad RP space onto these constrained 2 parameter spaces is complicated even in simple cases.
	* Chaos in the Dynamical System
	* Time Integrator Inaccuracy
	* Model Identifiability
	* Global Optimization
* Production models are simplified places which are easier to hunt down the many computational issues, and are simple enough to make it possible to understand the mechanisms 
* At the link provided here you can see our anlysis of the mechanisms of Bias for the Schaefer Model.


#Schnute 1985

* 

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

