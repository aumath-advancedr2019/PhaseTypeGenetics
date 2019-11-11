# PhaseTypeGenetics

The aim of this package is to provide the user with some comprehensive tools to analyse discrete and continuous phase-type distributions. As recent research confirms the applicability of phase-type theory in population genetics, most of the examples provided in the package are based on coalescent theory. In addition, the package provides some function that are intended for the use in population genetics only.

All functions and applications are based on the theory developed in 

* [BN] Mogens Bladt and Bo Friis Nielsen (2017): 
  *Matrix-Exponential Distributions in Applied Probability*. 
  Probability Theory and Stochastic Modelling (Springer), Volume 81 
* and [HSB] Asger Hobolth, Arno Siri-Jégousse, Mogens Bladt (2019): 
  *Phase-type distributions in population genetics*. 
  Theoretical Population Biology, 127, pp. 16-32.

To learn more about phase-type distributions and their applications in population genetics consider `vignette("PhaseTypeGenetics")`.

### Installation

```{r, eval=FALSE}
## Installing version from GitHub

devtools::install_github("aumath-advancedr2019/phastypdist/tree/master/PhaseTypeGenetics")

```

### Usage

The package PhaseTypeGenetics can be used to analyse phase-type distributions in general, but it also includes some functions that are aimed at the application in population genetics. All of the provided functions can only be used on objects of class `discphasetype` and `contphasetype`, hence it is necessary to use the constructors provided with the package.
