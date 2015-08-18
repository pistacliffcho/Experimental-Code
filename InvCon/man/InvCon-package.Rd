\name{InvCon-package}
\alias{InvCon-package}
\alias{InvCon}
\docType{package}
\title{
Computing the Inverse Convex Estimator
}
\description{
Fits the inverse convex estimator with inverse.convex(). Estimated densities, probabilities and 
quantiles can be found for a given fit via dIC, pIC and qIC. 
}
\details{
\tabular{ll}{
Package: \tab InvCon\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-26\cr
License: \tab NA  \cr
}
Fits an estimate such it is restrained to being inverse convex, i.e. f(x) = 1/g(x) where g(x) is convex.
This fit is actually a local maximum, not a global: the likelihood function is unbounded. Much like
the gaussian mixture, although the likelihood is unbounded around a degenerate distribution,
a very useful estimator can be obtained from a local non-degenerate distribution. Unlike the gaussian
mixture, this local maximum appears unique. Furthermore, the domain of convergence to degenerate
mode (emprically) shrinks very quickly, such that this algorithm avoids that mode with very high probability 
as long there is at least 30 unique points in the dataset. 
}
\author{
Clifford Anderson-Bergman

Maintainer: Clifford Anderson-Bergman <cianders@uci.edu>
}
\references{
Semi- and Non-parametric Methods for Interval Censored
  Data with Shape Constraints}
Inverse Convex, Shape Constraints
\keyword{ Inverse Convex }

\examples{
  fit = inverse.convex(rnorm(500) )
  # Fits an inverse convex estimator to a sample
  
  plot(fit)    						
  # Plots the estimated CDF
  
  qIC(0.5, fit)  						
  # Estimates the median
}
