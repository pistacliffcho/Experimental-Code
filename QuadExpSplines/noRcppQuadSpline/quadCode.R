quadExp <- function(x, a, b, c){
	exp(x^2 * a + x * b + c)
}

quasiCloseForm <- function(x1, x2, a, b, c){
	sigma = 1/sqrt(-2 * a)
	mu = b * sigma^2
	r = c + mu^2/(2 * sigma^2)
	exp(r) * sqrt(2 * pi * sigma^2) * (pnorm(x2, mu, sigma) - pnorm(x1, mu, sigma) )
}