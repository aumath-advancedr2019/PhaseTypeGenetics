# phastypdist
A package containing density, distribution and quantile functions for phase-type distributions and various ways of simulating from these distributions.

Project manager : Jette   

Documentation executive : Jette

Quality executive : Peter


The package will provide some tools for analysing population genetics. More precisely, it will be based on the article Phase-type distribution in population genetics by Hobolth et al from 2019 and hence gather our knowledge gained during the identically named course held by Asger Hobolth in autumn 2019 at Aarhus University. 
We will use phase-type distributions (discrete and continuous) to describe models for coalescent (with and without recombination). This implies that the package will provide the following:

- Properties for a general phase-type distribution (density, distribution function, simulation, Laplace transform, quantiles, mean and variance)
- Summary() for phase-type distributions
- Discretization of continuous phase-type distribution.
- The probability function of the number of segregating sites in a sample of n sequences.
- The reward-transformed distribution
- The block-counting process ? 
- Example : Different summary statistics, their distribution, mean and variance 
- Estimation
- Sums and minimum/maximum of phase-type distributions 

To this day, there is no R package that uses phase-type distributions to analyse coalescent models. We want to make it possible for the user to use this new approach. We want to provide some tools for analysing population genetics by the aid of matrix manipulations, and hence make analysis faster. 
