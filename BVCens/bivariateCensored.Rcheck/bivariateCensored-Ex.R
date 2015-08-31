pkgname <- "bivariateCensored"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bivariateCensored')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bivariateNPMLE")
### * bivariateNPMLE

flush(stderr()); flush(stdout())

### Name: bivariateNPMLE
### Title: Computes the MLE for Bivariate Interval Censored Data
### Aliases: bivariateNPMLE

### ** Examples

  testData <- simBVCen()
  #simulate bivariate interval censored data
  
  bvcensFit <- bivariateNPMLE(testData)
  #Finds the MLE
  
  bvcensFit



cleanEx()
nameEx("optCliq")
### * optCliq

flush(stderr()); flush(stdout())

### Name: optCliq
### Title: Computes the MLE for a Binary Mixture Model
### Aliases: optCliq

### ** Examples

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



cleanEx()
nameEx("simBVCen")
### * simBVCen

flush(stderr()); flush(stdout())

### Name: simBVCen
### Title: Simulates Bivariate Interval Censored Data
### Aliases: simBVCen

### ** Examples

  testData <- simBVCen()
  #simulate bivariate interval censored data
  
 bvcenFit <- bivariateNPMLE(testData)
 
 bvcenFit



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
