\name{optCliq}
\alias{optCliq}
\docType{package}
\title{
Computes the MLE for a Binary Mixture Model
}
\description{
Computes the MLE for a Binary Mixture Model
}
\usage{
optCliq(cliqMat, tol = 10^-10, 
        inner_loops = 100, outer_loops = 20)
}
\arguments{
  \item{cliqMat}{n x m clique matrix. n = number of observations, m = number of components}
  \item{tol}{numerical tolerance}
  \item{inner_loops}{number of inner loops used}
  \item{outer_loops}{number of outer loops used}
}

\examples{
testData <- simBVCen()
#simulate bivariate interval censored data

cliqMat <- MLEcens::reduc(testData, cm = TRUE)$cm
#computes the cliqMat associated with data

cliqMat <- t(cliqMat)
#reduc returns an m x n matrix, so
#needs to be transposed for compatibility with optCliq

cliqFit <- optCliq(cliqMat)
#optimizes the component weights for clique matrix

cliqFit 
}
