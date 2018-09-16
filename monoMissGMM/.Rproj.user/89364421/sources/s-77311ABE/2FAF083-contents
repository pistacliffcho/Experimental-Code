# library(monoMissGMM)

## TESTING FORWARD SOLVE

nCol = 4
X = matrix(rnorm(nCol^2 * 2), ncol = nCol)
S = cov(X)
U = t(chol(S) )

res1 = forwardsolve(U, X[1,])
res2 = forwardSolve(U, X[1,], r = 0)
identical(res1, res2)



### TESTING OUTER PRODUCT
X = matrix( rnorm(9), nrow = 3) 
z = rnorm(3)

manual_X = X + z %*% t(z) * 0.5
add_weighted_outer_prod(z, X, 0.5)
identical(X, manual_X)
