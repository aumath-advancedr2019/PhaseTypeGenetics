---
title: "Requirements specification"
author: "Jette Steinbach and Peter Krarup"
date: "4 10 2019"
output: html_document
---

## The user interface

The package will include the following functions:

* dphasetype : the density function for a phase-type distribution with initial distribution initDist and subintensity matrix Tmat
* pphasetype : the distribution function for a phase-type distribution with initial distribution initDist and subintensity matrix Tmat
* rphasetype: simulating from a phase type distribution with initial distribution initDist and subintensity matrix Tmat
* qphasetype : The quantile function 
* MomentsPhaseType: the Laplace transform for a phase-type distribution with initial distribution initDist and subintensity matrix Tmat
* mean : Returns the mean for a a phase-type distribution with initial distribution initDist and subintensity matrix Tmat (mean() is a generic function)
* variance: : Returns the variance for a phase-type distribution with initial distribution initDist and subintensity matrix Tmat (is not a genric function)
* summary : prints the initial distribution and the subintensity matrix 
* RewTransDistribution : A function that computes the reward-transformed distribution for a phase-type distribution with original initial distribution initDist and subintensity matrix Tmat. Here, the reward vector has to be non-negative.
* DistSegregatingSites : The distribution of the number of segregating sites based on a phase-type distribution with initial distribution initDist and subintensity matrix Tmat. Furthermore, the mutation rate theta is needed.
* BlockCountProcess : A function to find the state space and the corresponding rate matrix for the block counting process for a given number n of samples in the standard coalescent ? 
* EstimationPhaseType : estimate the inital distribution and subintensity matrix from a phase-type distribution.
* DiscretPhaseType : Discretization of continuous phase-type distribution.
* sum : Sum of phase-type distributions
* min : Minimum of phase-type distributions
* max : Maximum of phase-type distributions

We will define two classes dphasetype (discrete) and cphasetype (continuous), in order to make it possible for the user to use the standard functions mean(), variance(), sum(), min() and max().


Furthermore, we will include a dataset which is composed of an initial distribution, a subintensity matrix and several reward vectors to be able to show examples of applications. And a dataset containing a Markov process in order to estimate parameters. 


