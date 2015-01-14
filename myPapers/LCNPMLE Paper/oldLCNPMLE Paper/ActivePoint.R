x = 1:7
betas = c(1, 2, 2, 2, 2, 2, 1)
plot(x, betas, ylab = expression(beta), type = "b", ylim = c(0.5, 3.5) )
prm_chg_ht(4, c(1,2,6,7), 1, x, betas)
lines(x, betas, lty = 2, type = "b")
points(x, betas)
legend("topright", c( expression(beta^{(t)}), expression(beta^{(t+1)}) ), lty = c(1,2) )