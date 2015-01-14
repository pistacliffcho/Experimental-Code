source("/Users/cliffordanderson-bergman/Desktop/LCNPMLE Paper/TumorData.R")

lc.ce.fit <- LC.NPMLE(l = CE.data[,1], u = CE.data[,2], display = FALSE, return.output = TRUE)
lc.ge.fit <- LC.NPMLE(l = GE.data[,1], u = GE.data[,2], display = FALSE, return.output = TRUE)

uc.ce.fit <- EMICM(CE.data)
uc.ge.fit <- EMICM(GE.data)

lc.ce.cdf.inf <- d2cdf(lc.ce.fit[[1]], lc.ce.fit[[2]])
lc.ge.cdf.inf <- d2cdf(lc.ge.fit[[1]], lc.ge.fit[[2]])

lc.ce.x <- lc.ce.cdf.inf[[1]]
lc.ce.cdf <- lc.ce.cdf.inf[[2]]
lc.ge.x <- lc.ge.cdf.inf[[1]]
lc.ge.cdf <- lc.ge.cdf.inf[[2]]

xlim <- c(min(c(lc.ce.x, lc.ge.x)), max(c(lc.ge.x, lc.ce.x) ) ) 

par(mfrow = c(1,2) ) 

plot(lc.ce.x, 1-lc.ce.cdf, xlim = xlim, xlab = "Time", ylab = "S(t)", type = "l", lty = 2, ylim = c(0,1),  main = "LC NPMLE")
lines(lc.ge.x, 1-lc.ge.cdf,lty = 3) 

lines(c(0, 2000), c(0.5, 0.5) )

legend("topright", legend = c("Conventional", "Germ Free"), lty = c(2, 3) )

plot(0, 0, xlim = xlim, xlab = "Time", ylab = "S(t)", type = "n", ylim = c(0,1), main = "Unconstrained NPMLE")
lines.UC(uc.ce.fit)
lines.UC(uc.ge.fit, 3)

lines(c(0, 2000), c(0.5, 0.5) )
