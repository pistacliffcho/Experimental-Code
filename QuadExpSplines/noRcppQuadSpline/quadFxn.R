firstQuad <- function(x, y){
	c = y[1]
	x_end = x[3] - x[1]
	a = ((y[3] - y[1]) - (y[2]-y[1]) * x_end/x[2]) / (x_end*x_end - x_end*x[2])
	b = (y[2] - y[1] - a * x[2] * x[2]) / x[2];
	
	coefList <- list(a = a, b = b, c = c)
}


evalQuad = function(x, coefList){
	return(x*x * coefList$a + x * coefList$b + coefList$c)
}
coefList = list(a = -1, b = 2, c = 2)
evalPoints = c(-6, -2, 0, 2, 6)
evalQuad(evalPoints, coefList)
x1s = c(-2, -1, 0)
ys = c(-6, -2, 0)

fit <- firstQuad(x1s - x1s[1], ys)

xEval = -20:0/10

plot(x1s, ys)
lines(xEval, evalQuad(xEval - x1s[1], fit))
fit

fullXs = c(-6, -2, 0, 2, 6)
fullYs = evalQuad(evalPoints, coefList)

plot(fullXs, fullYs)
fit1 = list(a = -1, b = 14, c = -46)
fit2 = fit1
fit3 = list(a = -1, b = 6, c = -6)
fit4 = list(a = -1, b = 2, c = 2)
fit5 = list(a = -1, b = -2,c = 2)
fit6 = list(a = -1, b = -10, c = 22)

x1 <- (fullXs[1]*10-10):(fullXs[1]*10)/10
x2  <- (fullXs[1]*10):(fullXs[2]*10)/10
x3  <- (fullXs[2]*10):(fullXs[3]*10)/10
x4  <- (fullXs[3]*10):(fullXs[4]*10)/10
x5  <- (fullXs[4]*10):(fullXs[5]*10)/10
x6  <- (fullXs[5]*10):(fullXs[5]*10 +10)/10

lines(x1, evalQuad(x1 - fullXs[1] , fit1))
lines(x2, evalQuad(x2 - fullXs[1] , fit2))
lines(x3, evalQuad(x3 - fullXs[2], fit3 ))
lines(x4, evalQuad(x4 - fullXs[3], fit4 ))
lines(x5, evalQuad(x5 - fullXs[4], fit5 ))
lines(x6, evalQuad(x6 - fullXs[5], fit6 ))

allxs = c(x1, x2, x3, x4, x5, x6)
lines(allxs, evalQuad(allxs, coefList), col = 'blue')