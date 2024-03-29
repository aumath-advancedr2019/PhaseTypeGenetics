% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DocumentationBlockCountingProcess.R
\name{BlockCountProcess}
\alias{BlockCountProcess}
\title{Block counting process}
\usage{
BlockCountProcess(n)
}
\arguments{
\item{n}{the sample size (>=1)}
}
\value{
The function returns a list containing the sub-intensity
rate matrix \code{Rate_Mat} and the state space matrix \code{StateSpace_Mat}.
In the latter, each row corresponds to a state and each state is a
\eqn{(n-1)}-dimensional row vector.
}
\description{
Computing the state space and the corresponding rate matrix for
the block counting process for a given sample size \eqn{n} in
the standard coalescent model.
}
\details{
For a given sample size \eqn{n}, one can have one or more possible
coalescent trees. Each coalescent event in these trees correspond
to a state of the block counting process. Furthermore, each state is
represented by a \eqn{(n-1)}- dimensional row vector, where each entry \eqn{i}
corresponds to the number of branches giving rise to \eqn{i}
descendants. Hence, state 1 is always a vector of the form \eqn{(n,0,0,...,0)},
and state 2 is always given by the vector \eqn{(n-2,1,0,...,0)} eqn{(n >= 3)}.
}
\examples{
a <- BlockCountProcess(4)
a$Rate_Mat
a$StateSpace_Mat

}
