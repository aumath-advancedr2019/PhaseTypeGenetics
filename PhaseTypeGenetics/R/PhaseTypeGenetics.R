#' PhaseTypeGenetics: A package providing several properties of the phase-type
#' distribution with applications in population genetics
#'
#' The PhaseTypeGenetics package provides several functions to compute
#' properties of the discrete and continous phase-type distribution.
#' These functions are
#' \itemize{
#' \item contphasetype
#' \item discphasetype
#' \item dphasetype
#' \item maxima
#' \item mean
#' \item minima
#' \item moments
#' \item pphasetype
#' \item qphasetype
#' \item rphasetype
#' \item sum
#' \item summary
#' \item var
#' }
#' Furthermore, the package includes some functions that are useful in
#' connection to population genetics. These functions are
#' \itemize{
#' \item BlockCountProcess
#' \item discretization
#' \item dSegregatingSites
#' \item RewTransDistribution
#' \item SiteFrequencies,
#' }
#' The function \code{BlockCountProcess} computes the state space
#' and the corresponding rate matrix for
#' the block counting process for a given sample size \eqn{n} in
#' the standard coalescent model. The function \code{discretization}
#' discretizes a continuous phase-type distribution, and can for example
#' be used to obtain the distribution of the number of segregating sites,
#' when the mutation rate is known and equal to \eqn{\lambda}.
#' \code{dSegregatingSites} computes the distribution of the number
#' of segregating sites for a given sample size and a given mutation rate.
#' \code{RewTransDistribution} is an implementation of the reward
#' transformation and can be used to obtain the distribution of the
#' total branch length. Finally, the function \code{Site Frequencies}
#' computes the distribution of the site frequencies, the number of
#' segregating sites and the tail statistics using phase-type theory.
#'
#' @section PhaseTypeGenetics general functions:
#' \itemize{
#' \item \code{\link{contphasetype}}
#' \item \code{\link{discphasetype}}
#' \item \code{\link{dphasetype}}
#' \item \code{\link{maxima}}
#' \item \code{\link[=mean.discphasetype]{mean}}
#' \item \code{\link{minima}}
#' \item \code{\link{moments}}
#' \item \code{\link{pphasetype}}
#' \item \code{\link{qphasetype}}
#' \item \code{\link{rphasetype}}
#' \item \code{\link[=sum.discphasetype]{sum}}
#' \item \code{\link[=summary.discphasetype]{summary}}
#' \item \code{\link{var}}.
#' }
#'
#' @section PhaseTypeGenetics functions for applications in population genetics:
#' \itemize{
#' \item \code{\link{BlockCountProcess}}
#' \item \code{\link{discretization}}
#' \item \code{\link{dSegregatingSites}}
#' \item \code{\link{RewTransDistribution}}
#' \item \code{\link{SiteFrequencies}}.
#' }
#'
#' @docType package
#' @name PhaseTypeGenetics
#'
NULL

