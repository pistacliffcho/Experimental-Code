# library(monoMissGMM)

## CHECKING THE dmvnorm_one function. This isn't really used!!
nVars = 4
x = rnorm(nVars)
mu = rnorm(nVars)
X = matrix(rnorm(nVars^2 * 2), ncol = nVars)
S = cov(X)

# S = diag(nVars) * 2

chol_s = t(chol(S))
#diag(chol_s) = 1

monoMissGMM::dmvnorm_one(x, mu, chol_s, 0)
LaplacesDemon::dmvnc(x, mu, t(chol_s), log = T) 

# Testing dmvnorm_Omega (will be used)
split_res = monoMissGMM::dataPrep(X)
dmvnorm_Omega(split_res$data, mu, chol_s, split_res$Omega)
LaplacesDemon::dmvnc(X, mu, t(chol_s), log = T)

X[1,nVars] = NA
split_res = monoMissGMM::dataPrep(X)
dmvnorm_Omega(split_res$data, mu, chol_s, split_res$Omega)
LaplacesDemon::dmvnc(X[1,-nVars], mu[-nVars],
                     t(chol_s)[-nVars, -nVars], log = T)
LaplacesDemon::dmvnc(X[-1,], mu,
                     t(chol_s), log = T)
#...looks good