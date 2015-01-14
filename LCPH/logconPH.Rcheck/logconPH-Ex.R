pkgname <- "logconPH"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('logconPH')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dLC")
### * dLC

flush(stderr()); flush(stdout())

### Name: dLC
### Title: Density estimates from a log concave fit
### Aliases: dLC

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  dLC(0, fit)  						
  # Estimates the density at the true mode
  
  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  dLC(1, fit, covars = c(0,0))
  # Estimates the baseline density at t = 1



cleanEx()
nameEx("linesLC")
### * linesLC

flush(stderr()); flush(stdout())

### Name: linesLC
### Title: Draws Lines for Logconcave Fit
### Aliases: linesLC

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
  
  plotLC(fit)    					
  # Plots the estimated survival distribution
  
  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  plotLC(fit, covars = c(0,0), funtype = 'cdf')
  # Plots the estimated baseline cdf
  
  linesLC(fit, covars = c(1,1), funtype = 'cdf', col = 'red')
  # Plots the estimates cdf with covariates c(1,1)


cleanEx()
nameEx("logconPH-package")
### * logconPH-package

flush(stderr()); flush(stdout())

### Name: logconPH-package
### Title: Computing a Cox PH Model with a Log Concave Baseline
### Aliases: logconPH-package logconPH
### Keywords: Log Concave, CoxPH, Shape Constrained

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  qLC(0.5, fit)      				
  # Estimates the median

  simData <- sim_Censored(n = 400)
  # Simulates current status data
  
  fit = logconcave(simData)
  # Fits a log concave estimator to an interval censored sample

  pLC(0.5, fit)
  # Estimates the cdf at t = 0.5
  
  plotLC(fit,  'surv')
  # Plots the estimated survival function.
  # Options for second argument are 'pdf', 'cdf' and 'surv'

  simData <- simPH_Censored()
  # Simulates current status data from a Cox-PH model
  
  fit <- logconcave(times = simData$times, covariates = simData$x)
  # Fits a Cox-PH model with a logconcave baseline distribution
  
  plotLC(fit, covars = c(0,0) )
  # Plots the estimated baseline survival function
  
  linesLC(fit, covars = c(1,1), col = 'red')
  # Plots the estimated survival function with x1 = 1, x2 = 1



cleanEx()
nameEx("logconcave")
### * logconcave

flush(stderr()); flush(stdout())

### Name: logconcave
### Title: Cox PH model with Log Concave Baseline
### Aliases: logconcave

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  qLC(0.5, fit)      				
  # Estimates the median

  simData <- sim_Censored(n = 400)
  # Simulates current status data
  
  fit = logconcave(simData)
  # Fits a log concave estimator to an interval censored sample

  pLC(0.5, fit)
  # Estimates the cdf at t = 0.5
  
  plotLC(fit,  'surv')
  # Plots the estimated survival function.
  # Options for second argument are 'pdf', 'cdf' and 'surv'

  simData <- simPH_Censored()
  # Simulates current status data from a Cox-PH model
  
  fit <- logconcave(times = simData$times, covariates = simData$x)
  # Fits a Cox-PH model with a logconcave baseline distribution
  
  plotLC(fit, covars = c(0,0) )
  # Plots the estimated baseline survival function
  
  linesLC(fit, covars = c(1,1), col = 'red')
  # Plots the estimated survival function with x1 = 1, x2 = 1



cleanEx()
nameEx("pLC")
### * pLC

flush(stderr()); flush(stdout())

### Name: pLC
### Title: Probability estimates from a log concave fit
### Aliases: pLC

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  pLC(0, fit)  						
  # Estimates the cdf at the true mode

  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  pLC(1, fit, covars = c(0,0))
  # Estimates the baseline probability at t = 1



cleanEx()
nameEx("plotLC")
### * plotLC

flush(stderr()); flush(stdout())

### Name: plotLC
### Title: Plots Logconcave Fit
### Aliases: plotLC

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
  
  plotLC(fit)      				
  # Plots the estimated survival distribution
  
  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  plotLC(fit, covars = c(0,0), funtype = 'cdf')
  # Plots the estimated baseline cdf
  
  linesLC(fit, covars = c(1,1), funtype = 'cdf', col = 'red')
  # Plots the estimates cdf with covariates c(1,1)


cleanEx()
nameEx("qLC")
### * qLC

flush(stderr()); flush(stdout())

### Name: qLC
### Title: Quantiles estimates from a log concave fit
### Aliases: qLC

### ** Examples

  fit = logconcave(rnorm(500) )
  # Fits a log concave estimator to an uncensored sample
    
  qLC(0.5, fit)  						
  # Estimates the median

  simData <- simPH_Censored()
  # Simulates current status data from a CoxPH model
  
  fit <- logconcave(simData$times, simData$x)
  # Fits coxPH model
  
  qLC(0.5, fit, covars = c(0,0))
  # Estimates the baseline median



cleanEx()
nameEx("simCensored")
### * simCensored

flush(stderr()); flush(stdout())

### Name: simPH_Censored
### Title: Simulate current status data from Cox-PH model
### Aliases: simPH_Censored

### ** Examples

  simData <- simPH_Censored()
  # Simulates censored data from a Cox-PH model
  
  fit <- logconcave(times = simData$times, covariates = simData$x)
  # Fits a Cox-PH model with a logconcave baseline distribution



cleanEx()
nameEx("sim_Censored")
### * sim_Censored

flush(stderr()); flush(stdout())

### Name: sim_Censored
### Title: Simulate current status data from a beta(2,2) distribution
### Aliases: sim_Censored

### ** Examples

  simData <- sim_Censored()
  # Simulates current status data
  
  fit <- logconcave(simData)
  # Fits a log concave fit



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
