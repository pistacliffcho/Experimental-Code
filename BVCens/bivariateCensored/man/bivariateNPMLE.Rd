\name{bivariateNPMLE}
\alias{bivariateNPMLE}
\docType{package}
\title{
  Computes the MLE for Bivariate Interval Censored Data
}
\description{
  Computes the MLE for a Interval Censored Data
}
\usage{
  bivariateNPMLE(times, tol = 10^-10, 
          inner_loops = 100, outer_loops = 20)
}
\arguments{
  \item{times}{n x 4 matrix representing bivariate interval censored data. 
    The format is (l1, u1, l2, u2), i.e. the first column should be the lower limit of the first interval, 
    second column is upper limit of first interval, etc.}
  \item{tol}{numerical tolerance}
  \item{inner_loops}{number of inner loops used}
  \item{outer_loops}{number of outer loops used}
}

\examples{
  testData <- simBVCen()
  #simulate bivariate interval censored data
  
  bvcensFit <- bivariateNPMLE(testData)
  #Finds the MLE
  
  bvcensFit
}
