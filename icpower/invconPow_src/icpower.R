dyn.load('~/Desktop/Shape-Constraints/icpower/invconPow_src/logcon.so')

icpower = function(data){
	alpha = 1
	if(class(data) == "numeric"){
		x = sort(data)
		repVec <- as.numeric(table(x))
		x <- unique(x)
		inds = 1:length(x)
		if(max(x) == Inf | min(x) == -Inf)
			stop("Error: non finite values supplied")
		b = rep(1, length(x))
		output = .Call('uniVarICCens', inds, inds, x, b, repVec, as.integer( c(1, length(x) ) ), FALSE, as.numeric(alpha) )
		names(output) <- c("x", "beta", "llk")
		output[['dens']] <- (output$beta)^(-alpha)
		class(output) <- "ICPObject"
		return(output)
	}

if(is.matrix(data))
	int.Data = data
else
	stop("Data type not recognized")	
l = int.Data[,1]
u = int.Data[,2]
	

	
	if(max(l) <= min(u) )
		{
		x.use = c(max(l), min(u))
		density = c(1/(min(u) - max(l) ),1/(min(u) - max(l) ))
		lk.final = 0
		it = 0
		x = x.use
		l.dens = log(density)
		output <- list(x.use, density, lk.final, it, x, l.dens)
		cat("Note: Universal overlap leads to poorly performing estimator!\n")
		names(output) <- c("x", "dens", "lk", "it", "x.all", "l.dens")
	#	return(output)
		}	

	
	
	#NEED TO DECIDE ON WHAT TO DO ABOUT START VALUE FOR INVERSE CONVEX POWER...
	output = start.Inf(l, u)
	x = output[[1]]
	b = output[[2]]
	b = -b + 1
	repData <- rep.vec(int.Data)
	l <- repData[,1]
	u <- repData[,2]
	repVec = repData[,3]
	minX = if(x[1] == -Inf) 2 else 1
	maxX = if(x[length(x)] == Inf) length(x) - 1 else length(x) 
	actInds = c(minX, which(b == min(b) ) , maxX )
	actInds= unique(actInds)
	inds <- make.inds(x, l,u)
	output <- .Call('uniVarICCens', inds[,'l.inds'], inds[,'u.inds'], x, b, repVec, as.integer(actInds) , TRUE, as.numeric(alpha) )
	names(output) <- c("x", "beta", "llk")
	output[['dens']] <- (output$beta)^(-alpha)
	class(output) <- "ICPObject"
	return(output)
}


rep.vec <- function(cs.data)
	{
	drop <- min(cs.data[,1]) == cs.data[,1] & max(cs.data[,2]) == cs.data[,2]
	cs.data <- cs.data[!drop,]
	l <- cs.data[,1]
	u <- cs.data[,2]
	n <- length(l)
	min.l <- min(l)
	max.u <- max(u)
	left <- u[l == min.l]
	right <- l[u == max.u]	
	left <- sort(left)
	right <- sort(right)
	if(length(left) + length(right) != n)
		{
		rep.vec.output <- makeNonCS_RepVec(l, u)
		return(cbind(rep.vec.output[[1]], rep.vec.output[[2]], rep.vec.output[[3]], row.names = NULL) ) 
		}
	rep.lft <- table(left)
	val.lft <- unique(left)
	rep.rgt <- table(right)
	val.rgt <- unique(right)
	l <- c( rep(min.l, length(val.lft)), val.rgt )
	u <- c( val.lft, rep(max.u, length(val.rgt) ) )
	rep.vec <- c(rep.lft, rep.rgt) 
	names(rep.vec) <- NULL
		return(cbind(l, u, rep.vec,  row.names = NULL) )
	}		
	
make.inds <- function(x, l, u)
	{
	u.inds <- NA
	l.inds <- NA
	for(i in 1:length(l) ) 
		{
		u.inds[i] <- which(x == u[i])
		l.inds[i] <- which(x == l[i])
		}
	return(cbind(u.inds, l.inds) ) 
	}
	
	
start.Inf <- function(l, u)
	{
	x.pos <- sort( unique( c(l, u) ) )
	
	k = length(x.pos)
	
	x.mid = (x.pos[-1] + x.pos[-k])/2
	if(x.pos[1] == -Inf)
		x.mid[1] = x.pos[2] - 1
	if(x.pos[k] == Inf)
		x.mid[k-1] = x.pos[k-1] + 1
	x <- sort( c(x.pos, x.mid) )
	k <- length(x)
	
	mid.ind = round(k/2)
	x.mid = x[mid.ind]
	
	
	betas = -abs(x - x.mid)/sd(x[abs(x) < Inf])
	
	return(list(x, betas) )
	}
