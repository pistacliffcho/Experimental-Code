library(quadprog)
library(Rcpp)

#sourceCpp("~/Desktop/MySoftware/InverseConvex/InverseConvex1.0.cpp")


# High Level Functions:

#		inverse.convex(x, max.it = 500, tol = 10^-4, plot = FALSE, display = FALSE)			
#		Fits the inverse convex NPMLE
#		Arguments:
#		x			data set to be fit with inverse convex estimator
#		max.it		maximum iterations of algrorithm
#		tol			tolerance required for convergence
#		plot		should the estimated density be plotted? 
#		display		should information about the convergence of the algorithm be displayed?
#		Value:
#		Returns an "ICfit" object, which contains "x", the location of the active points and "dens",
#		the estimated density at the active points. These are the sufficient statistics for describing
#		an inverse convex fit

#		Example:
#		fit = inverse.convex(rnorm(500) )  		# Fits an inverse convex estimator to a sample
#		plot(fit)								# Plots the estimated CDF
#		qIC(0.5, fit)							# Estimates the median


#		plot.ICfit(ICfit, type = "cdf", ...)	
#		Plots the estimated density function for an inverse convex fit
#		Arguments:		
#		ICfit			an ICfit object returned from the inverse.convex function
#		type			type of function plotted. Current choices: "cdf", "pdf"
#		...				additional arguments to be passed to the plot function	

#		pIC(q, ICfit)		
#		Returns estimated P(X > q) for inverse convex fit
#		Arguments:
#		q		quantile
#		ICfit	ICfit object
#		Values:
#		returns a vector of estimated probabilities	

#		dIC(x, ICfit)		
#		Returns estimated density for inverse convex fit
#		Arguments:
#		x		value at which density is calculated
#		ICfit	ICfit object
#		Values:
#		returns a vector of estimated densities	

#		qIC(p, ICfit)		
#		Returns estimated quantiles for inverse convex fit
#		Arguments:
#		p		probabilities
#		ICfit	ICfit object
#		Values:
#		returns a vector of estimated quantiles	

#'		Computes the Inverse Convex Estimator
#'
#'		
inverse.convex = function(x, max.it = 1000, tol = 10^-4, plot = FALSE, display = FALSE, low.lim = 10^-6, max.d = 10^7, start = "unif")
	{
	if(sum(is.na(x) ) > 0 )
		stop("NA's in data provided to inverse.convex")
	# Computes Inverse convex estimator
	sd.est = sd(x)
	x = sort(x)
	w = as.vector(table(x) )
	x = unique(x) 
	k = length(x)
	betas = rep(5, k) 
	it = 0
	err = tol + 1
	lk.old = ic_llk(x, betas, w)
	x = x/sd.est

	if(start == "right")
		betas = betas + (max(x) - x)
	if(start == "left")
		betas = betas + (x - min(x))
		
	PROB = FALSE
	while(err > tol & it < max.it & PROB == FALSE)
		{
		it = it + 1
		update = ic.vdm(x, betas, w, low.lim = low.lim)
		err = update[[2]]
		if(err > 10^20)
			{
			cat("warning: err > 10^20! Most likely approaching degenerate distribution\n")
			PROB = TRUE
			} else{
			betas = update[[1]]
			betas = ic.acts.opt(x, betas, w, low.lim = low.lim)
			betas = ic.ICM(x, betas, w, low.lim = low.lim)
			betas = betas * ic_mass(x, betas)
			if(min(betas) < 10^-10 & PROB == FALSE)
				{
				PROB = TRUE
				cat("warning: max estimated density > max.d! Most likely approaching degernerate distribution\n")
				}
			}
		}
	x = x * sd.est
	lk.new = ic_llk(x, betas, w)
	betas = betas * ic_mass(x, betas) 	
	if(display == TRUE)
		{
		cat("err = ", err, " iterations = ", it, "\n")
		cat("lk.old = ", lk.old, "lk.new = ", lk.new, "\n")
		}
	mass = ic_mass(x, betas)
	betas = betas * mass
	
	if(plot == TRUE)
		{
		plot(x, 1/betas, type = "l", xlab = "x", ylab = "Estimated Density")
		act.inds = ic.active(x, betas)
		points(x[act.inds], 1/betas[act.inds])
		}
		
	active = ic.active(x, betas)			
	x = x[active]
	betas = betas[active]
	mass = ic_mass(x, betas)
	dens = 1/betas
	dens = dens/mass
	if(err > tol)
		PROB = TRUE
	output = list(x, dens, ic_llk(x, betas, w), !PROB ) 
	names(output) = c("x", "density", "lk", "conv")
	class(output) = "ICfit"
	return(output )
	}
	
plot.ICfit <- function(ICfit, type = "cdf", ...)
	{
	#plots ICfit object
	grid = min(ICfit$x) + 1:999/1000 * (max(ICfit$x) - min(ICfit$x))
	y = NA
	if(type == "pdf")
		{
		for(i in 1:length(grid) ) 
			y[i] = d_ic(grid[i], ICfit$x, ICfit$dens)
		plot(grid, y, xlab = "x", ylab = "Estimated Density", type = "l", ...)
		}
	if(type == "cdf")
		{
		for(i in 1:length(grid) ) 
			y[i] = p_ic(grid[i], ICfit$x, ICfit$dens)
		plot(grid, y, xlab = "x", ylab = "Estimated CDF", type = "l", ...)
		}
	}	
	
	
find.limits = function(ind, act.inds, x, betas, low.lim = 10^-4, max.fact = 5)
	{
	# Finds limits for active points dictated by shape constraint	
	l.lim = low.lim
	h.lim = betas[ind] * max.fact
	act.inds = sort(unique( c(ind, act.inds) ) )
	a.ind = which(ind == act.inds)
	ak = length(act.inds)
	if(ak == 2)
		return(c(l.lim, h.lim) ) 
	if(a.ind > 1 & a.ind < ak)
		{
		dx.i1 = x[ind] - x[act.inds[a.ind-1] ] 
		dx.i = x[act.inds[a.ind+1]] - x[ind]
		h.lim = min( c(h.lim, (betas[act.inds[a.ind+1]]/dx.i + betas[act.inds[a.ind-1] ]/dx.i1) / (1/dx.i + 1/dx.i1) ) )
		}
	if(a.ind > 2)
		{
		dx.i1 = x[ind] - x[act.inds[a.ind-1] ] 
		dx.i2 = x[act.inds[a.ind-1] ] - x[act.inds[a.ind-2] ] 
		db.i2 = betas[act.inds[a.ind-1] ] - betas[act.inds[a.ind-2] ] 
		l.lim = max( c(l.lim, betas[act.inds[a.ind-1]] + dx.i1 * db.i2/dx.i2) )
		}
	if(a.ind < (ak - 1) )
		{
		dx.i = x[act.inds[a.ind+1]] - x[ind]
		dx.i1 = x[act.inds[a.ind+2]] - x[act.inds[a.ind+1]]
		db.i1 = betas[act.inds[a.ind+2]] - betas[act.inds[a.ind+1]]
		l.lim = max(c(l.lim, betas[act.inds[a.ind+1] ] - dx.i * db.i1 / dx.i1 ) )
		}
	return(c(l.lim, h.lim) ) 
	}
	
ic.active <- function(x, betas, slack = 10^-10)
	{
	# Finds active points of current fit
	active = rep(FALSE, length(x) ) 
	ic_active(active, x, betas, slack)
	act.inds = which(active)
	ak = length(act.inds)
	if(ak > 2)
		{
		drop = rep(FALSE, ak-1)
		
		slopes = (betas[act.inds[-1]] - betas[act.inds[-ak]])/(x[act.inds[-1]] - x[act.inds[-ak]])
		for(i in 1:(ak-2) ) 
			{
			drop[i] = slopes[i+1] - slack < slopes[i]
			}
		if(sum(drop) > 0)
			{
			act.inds = act.inds[-(which(drop) + 1)]
			}
		}
	return(act.inds)
	}

ic.uni.opt <- function(ind, x, betas, weight, h = 10^-4, low.lim)
	{
	# This function both provides univariate optimization of a given active point
	low.lim = betas[ind]/1.5
	betas.use = betas
	betas.use[1] = betas.use[1] + 0
	lk.old <- ic_llk(x, betas.use, weight)
	ders <- c(0,0)
	act.inds = ic.active(x, betas.use)
	ic_act_ders(ind, ders, act.inds, x, betas.use, weight, h)
	d <- ders[1]	
	d2 <- ders[2]
	min.d <- -Inf
	max.l <- Inf
	max.r <- Inf
	limits = find.limits(ind, act.inds, x, betas.use, low.lim = low.lim)
	min.d = limits[1] - betas.use[ind]
	max.d = limits[2] - betas.use[ind]
	if(is.na(d2))
		{
		cat("warning: is.na(d2) in VEM step! Summary of betas...\n")
		print(summary(betas.use[betas.use > -Inf]) ) 
		cat("err = ", d, "\n")
		cat("lk = ", lk.old, "\n")
		cat("min(betas) = ", min(betas), "\n")
		}
	if(d2 < 0)
		{
		dlt <- -d/d2
		dlt <- max(min.d, dlt)
		dlt <- min(dlt, max.d)
		if(dlt == 0)
			{
			return(betas.use)
			}
		move_active(ind, dlt, act.inds, x, betas.use)
		lk.new <- ic_llk(x, betas.use, weight) 
		if(lk.new > lk.old)
			{
			return(betas.use)
			}
		for(i in 1:5)
			{
			dlt <- dlt/2
			move_active(ind, -dlt, act.inds, x, betas.use)
			lk.new <- ic_llk(x, betas.use, weight) 
			if(lk.new > lk.old)
				return(betas.use)
			}
		}
	dlt <- sign(d) * 4
	dlt <- max(min.d, dlt)
	dlt <- min(max.d, dlt)
	if(dlt == 0)
		{
		return(betas.use)
		}
	move_active(ind, dlt, act.inds, x, betas.use)
	it = 0
	dlt <- -dlt/2
	lk.new <- ic_llk(x, betas.use, weight) 
	while(it < 20 & lk.new < lk.old)
		{
		it <- it + 1
		move_active(ind, dlt, act.inds, x, betas.use)
		lk.new <- ic_llk(x, betas.use, weight) 
		if(lk.new > lk.old)
			return(betas.use)
		dlt = dlt/2
		}
	if(lk.new > lk.old)
		return(betas.use)	
		
	return(betas)
	}		
	
ic.vdm <- function(x, betas, weight, h = 10^-4, low.lim)
	{
	#This function finds the point corresponding to the highest KKT error and performs univariate optimization on it
	active = rep(FALSE, length(x) ) 
	ic_active(active, x, betas)	
	act.inds= which(active)

	ind.err = rep(0, length(x) )
	ic_aerr_vec(ind.err, act.inds, x, betas, weight)
	err.vec = rep(0, length(x) ) 
	act_d_vec(x, ind.err, act.inds, err.vec)
		
	err.vec2 = rep(0, length(x) ) 	
		
	err.vec[!active & err.vec > 0] = 0

	err = max(abs(err.vec) )
	
	if(!is.numeric(err) ) 
		{
		cat("Warning: algorithm error undefined. May be approaching degenerate distribution\n")
		err = Inf
		return(list(betas, err) ) 
		}
	
	err.ind = which(abs(err.vec) == err)[1]
	
	if(is.na(err.ind) ) 
		{
		cat("Warning: algorithm error undefined. May be approaching degenerate distribution\n")
		err = Inf
		return(list(betas, err) ) 
		}



	betas = ic.uni.opt(err.ind, x, betas, weight, h, low.lim)
	return(list(betas, err) ) 
	}

ic.acts.opt <- function(x, betas, weight, h = 10^-4, low.lim)
	{
	#Performs univariate optimization on all current active points
	org.act.inds = ic.active(x, betas)
	ak = length(org.act.inds)
	for(i in 1:ak)
		betas = ic.uni.opt(org.act.inds[i], x, betas, weight, low.lim = low.lim)
	return(betas)
	}


ic.cnst <- function(betas, x, act.inds)
	{
	# Finds the current constrant matrix and vector required for quadratic programming
	x <- x[act.inds]
	betas <- betas[act.inds]
	k <- length(betas)
	dx <- x[-1] - x[-k]
	con.vec <- rep(0, k) 
	for(i in 1:(k-2))
		con.vec[i] <- (betas[i+2] - betas[i+1])/dx[i+1] - (betas[i+1] - betas[i])/dx[i]
	con.vec[k-1] <- con.vec[k-2]
	con.vec[k] <- con.vec[k-1]
	
	con.vec = -con.vec 
	
	con.mat <- matrix(0, nrow = k, ncol = k)
	for(i in 1:(k-2) )
		{
		con.mat[i,i] <- 1/dx[i]
		con.mat[i, i+1] <- - 1/dx[i] - 1/dx[i+1]
		con.mat[i, i+2] <- 1/dx[i+1]
		}
	con.mat[k-1,] <- con.mat[k-2,]
	con.mat[k,] <- con.mat[k-1,]
	con.mat <- t(con.mat)

	return(list(con.mat, con.vec) ) 
	}	

	
ic.ICM <- function(x, betas, w, h = 10^-4, thresh = 1, low.lim)
	{
	# Performs ICM step, i.e. updates all active points simulatenously 
	act.inds <- ic.active(x, betas)
	a.k <- length(act.inds)
	if(a.k == 2)
		return(betas)
	d.vec <- NA
	d2.vec <- NA
	lk.old <- ic_llk(x, betas, w)
	i <- 0
	ders <- c(0,0)
	drop = rep(FALSE, a.k)
	while(i < a.k)
		{
		i <- i + 1
		a.ind <- act.inds[i]
		ic_act_ders(a.ind, ders, act.inds, x, betas, w)
		d.vec[i] <- ders[1]
		d2.vec[i] <- ders[2]
		if(ders[2] >= 0 )
			return(betas)
		ders <- c(0,0)
		}
	ind.err = rep(0, length(x) )
	ic_aerr_vec(ind.err, act.inds, x, betas, w)
	err.vec = rep(0, length(x) ) 
	act_d_vec(x, ind.err, act.inds, err.vec)
	d.vec = err.vec[act.inds]	
	con.inf <- ic.cnst(betas, x, act.inds)
	con.mat <- con.inf[[1]]
	con.vec <- con.inf[[2]]
	max.con <- max(abs(con.mat) ) * 10
	if(max.con > thresh)
		{
		con.mat <- con.mat/max.con
		con.vec <- con.vec/max.con
		}
		
	d2.mat <- matrix(0, nrow = a.k, ncol = a.k)
	diag(d2.mat) <- d2.vec
	con.vec <- con.vec - 10^-12
	
	
	delta <- rep(0, length(act.inds) )
	
	min.beta = min(betas)
	
	if(min(abs(con.vec) ) > 10^-12)
		{	
		delta <- try(solve.QP(Dmat = -d2.mat, dvec = -d.vec, Amat = -con.mat, bvec = con.vec)$solution, silent = TRUE)
		if(!is.vector(delta) ) 
			return(betas)
		delta <- -delta
		}
	
	for(i in 1:a.k)
		{
		if( -delta[i] > betas[act.inds[i]] )
			delta = delta * 3/4 * betas[act.inds[i]]/delta[i]
		}
		
	for(i in 1:a.k)
		{
		move_active(act.inds[i], delta[i], act.inds, x, betas)
		}
	hf.stp <- 0
	max.stp <- 10
	
	delta = -delta
	
	lk.new <- ic_llk(x, betas, w)
	if(is.na(lk.new) ) 
		lk.new = -Inf		
	while(lk.new < lk.old & hf.stp < max.stp)
		{
		hf.stp <- hf.stp + 1
		delta <- delta/2
		for(i in 1:a.k)
			{
			move_active(act.inds[i], delta[i], act.inds, x, betas)
			}
		lk.new <- ic_llk(x, betas, w)
	if(is.na(lk.new) ) 
		lk.new = -Inf		
		}
	if(hf.stp == max.stp)
		{
		move_active(act.inds[i], delta[i], act.inds, x, betas)
		}
	if(lk.new > lk.old)
		return(betas)
	return(betas)
	}	
	
pIC = function(q, ICfit)
	{
	#Estimated CDF
	p = NA
	for(i in 1:length(q) ) 
		p[i] = p_ic(q[i], ICfit$x, ICfit$dens)
	return(p)
	}	
	
dIC = function(x, ICfit)
	{
	#Estimated PDF
	d = NA
	for(i in 1:length(x))
		d[i] = d_ic(x[i], ICfit$x, ICfit$dens)
	return(d)
	}

qic.opt.fxn = function(q, p0, ICfit)
	{
	#Required for estimated quantile
	err = (p_ic(q, ICfit$x, ICfit$dens) - p0)^2
	return(err)
	}

qIC = function(p, ICfit)
	{
	#Estimated quantile
	q = NA
	inter = c(ICfit$x[1], ICfit$x[length(ICfit$x)])
	for(i in 1:length(p) ) 
		q[i] = optimize(qic.opt.fxn, inter, p0 = p[i], ICfit = ICfit)$minimum
	return(q)
	}