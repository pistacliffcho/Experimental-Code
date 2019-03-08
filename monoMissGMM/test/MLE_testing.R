library(monoMissGMM)

nRow = 1000
nCol = 4

dat = matrix(rnorm(nRow * nCol, mean = 1:4, sd = 1:4),
             ncol = nCol, byrow = T)
dat[1,(nCol-2):nCol] = NA
dat[2,nCol] = NA


info = monoMissGMM:::dataPrep(dat)

ws = runif(nRow)
sum_yw(info$data, info$Omega, 0, ws) / sum(ws[-(1:2)])
# The following should be equal to the 
# first elements of the above vector
weighted.mean(dat[-(1:2), 1], w = ws[-(1:2)])
weighted.mean(dat[-(1:2), 2], w = ws[-(1:2)])

monoMissGMM:::compute_y_r(info$data, info$Omega, ws)
# First element should be equal to the vector above
# Below should be equal to the last element
weighted.mean(dat[,1], w = ws)


q = ncol(dat)
ss_start = diag(q) * 0
mle = monoMissGMM::computeMLE(info$dat, ws, info$Omega, ss_start)
mle$mu_hat
mle$chol_hat
mle_vec = dmvnorm_Omega(info$data, mle$mu_hat, mle$chol_hat, info$Omega)
eps = 0.001
up_vec = dmvnorm_Omega(info$data, mle$mu_hat + eps, mle$chol_hat, info$Omega)
dn_vec = dmvnorm_Omega(info$data, mle$mu_hat - eps, mle$chol_hat, info$Omega)

sum(mle_vec * ws)
sum(up_vec * ws)
sum(dn_vec * ws)
mle$mu_hat



library(gmmMonoMiss)

data2 = dat[,nCol:1]
split_info = gmmMonoMiss::makeDataList(data2, cum = F)
n_missings = get_split_patterns(split_info$split_data)
split_w = split(ws, split_info$missing_vec)
mle2 = gmmMonoMiss::compute_mle(split_info$split_data, 
                                split_w, n_missings)

mu_org = mle2$mu[nCol:1]
chol_org = mle2$chol[nCol:1, nCol:1]
