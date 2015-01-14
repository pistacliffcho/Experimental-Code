n = 250
data <- rexp(n)
fit <- logconcave(data)
par(mfrow = c(1,2))
x <- min(data) + 0:100/100 * (max(data) - min(data))
yLims = c(0,5)
plot(x, dLC(x, fit)/(1 - pLC(x, fit)), xlab = 't', ylab = 'Estimated h(t)', main = 'Log Concave Estimated Hazard', type = 'l', ylim = yLims)
lines(x, dexp(x)/(1-pexp(x)), col = 'red')
augD <- (n-1)/n
augS <- .5/n
plot(x, dLC(x, fit) * augD/(1 - pLC(x, fit) * augD + augS), xlab = 't', ylab = 'Estimated h(t)', main = 'Augmented\nLog Concave Estimated Hazard', type = 'l', ylim = yLims)
lines(x, dexp(x)/(1-pexp(x)), col = 'red')
legend('topright', legend = c('Estimate', 'True Hazard'), lwd = c(1,1), col = c('black', 'red'))