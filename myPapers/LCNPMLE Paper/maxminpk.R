require(graphics)
matplot(0, 0, xlab = "x", ylab = expression(paste(phi, "(x)") ), type = "n", xlim = c(1, 7), ylim = c(0,4), xaxt = "n") 

lines( c(1,2), c(0,1) )
lines(c(6,7), c(1,0) ) 

lines( c(2, 6), c(1,1), lty = 3) 
lines( c(2, 4), c(1, 3), lty = 2)
lines(c(4,6), c(3,1), lty = 2)

lgd <-  expression( "fixed", paste("max ", p[i]),  paste("min ", p[i])  ) 
legend("topright", legend = lgd, lty = c(1, 2, 3) )

axis(side = 1, at = c(2, 4, 6), labels = expression(l[i], m[i], r[i]) )