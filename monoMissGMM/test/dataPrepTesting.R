# library(monoMissGMM)

nRow = 5
nCol = 4

dat = matrix(rnorm(nRow * nCol), ncol = nCol)
dat[1,(nCol-2):nCol] = NA
dat[2,nCol] = NA


res = monoMissGMM:::dataPrep(dat)
res
# Should have data matrix, and Omega list
# Omega list should have 0-based indices for each missingness pattern