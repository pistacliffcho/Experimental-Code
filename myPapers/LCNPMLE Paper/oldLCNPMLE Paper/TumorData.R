CE.no.t <- c(45, 198, 215, 217, 257, 262, 266, 371, 431, 447, 454, 459, 475, 479, 484, 500, 502, 503, 505, 508, 516, 531, 541, 553, 556, 570, 572, 575, 577, 585, 588, 594, 600, 601, 608, 614, 616, 632, 638, 642, 642, 642, 644, 644, 647, 647, 653, 659, 660, 662, 663, 667, 667, 673, 673, 677, 689, 693, 718, 720, 721, 728, 760, 773, 777, 815, 886)

CE.t <- c(381, 477, 485, 515, 539, 563, 565, 582, 603, 616, 624, 650, 651, 656, 659, 672, 679, 698, 702, 709, 723, 731, 775, 779, 795, 811, 839)

n.ce.nt <- length(CE.no.t)
n.ce.t <- length(CE.t)

ce.u <- c(CE.t, rep(1100, n.ce.nt) )
ce.l <- c(rep(0, n.ce.t), CE.no.t )

GE.t <- c(546, 609, 692, 692, 710, 752, 773, 781, 782, 789, 808, 810, 814, 842, 846, 851, 871, 873, 876, 888, 888, 890, 894, 896, 911, 913, 914, 914, 916, 921, 921, 926, 936, 945, 1008)

GE.nt <- c(412, 524, 647, 648, 695, 785, 814, 817, 851, 880, 913, 942, 986)

n.ge.nt <- length(GE.nt)
n.ge.t <- length(GE.t)

ge.u <- c(GE.t, rep(1100, n.ge.nt) )
ge.l <- c(rep(0, n.ge.t), GE.nt )

GE.data <- cbind(ge.l, ge.u)
CE.data <- cbind(ce.l, ce.u) 

Ho.samp <- function(CE.data, GE.data)
	{
	n.ce <- length(CE.data[,1])
	n.ge <- length(GE.data[,1])
	
	comb.data <- rbind(CE.data, GE.data)
	samp <- sample(1:(n.ce + n.ge), n.ge)
	
	ge.samp <- comb.data[samp,]
	ce.samp <- comb.data[-samp, ]
	return(list(ge.samp, ce.samp) ) 
	}
	
per.vals <- function(CE.data, GE.data, MC =	10)
	{
	uc.ge.med <- NA
	lc.ge.med <- NA
	uc.ce.med <- NA
	lc.ce.med <- NA
	for(i in 1:MC)
		{
		perm.samp <- Ho.samp(CE.data, GE.data)
		ge.samp <- perm.samp[[1]]
		ce.samp <- perm.samp[[2]]
		uc.ge <- EMICM(ge.samp)
		uc.ce <- EMICM(ce.samp)
		uc.ge.med[i] <- uc.med(uc.ge)
		uc.ce.med[i] <- uc.med(uc.ce)
		lc.ge <- LC.NPMLE(l = ge.samp[,1], u = ge.samp[,2], return.out = TRUE, display = FALSE)
		lc.ce <- LC.NPMLE(l = ce.samp[,1], u = ce.samp[,2], return.out = TRUE, display = FALSE)
		lc.ge.med[i] <- lc.med(lc.ge)
		lc.ce.med[i] <- lc.med(lc.ce)
		}
	return(list(uc.ge.med, uc.ce.med, lc.ge.med, lc.ce.med) )
	}
	
sd.pe.samp <- function(pe.vals)
	{
	cat("SD of uc.ge.med = ", sd(pe.vals[[1]]),"\nSD of lc.ge.med = ", sd(pe.vals[[3]]), "\nSD of uc.ce.med = ", sd(pe.vals[[2]]), "\nSD of lc.ce.med = ", sd(pe.vals[[4]]) )
	
	cat("\nSD of uc.ge.med - uc.ce.med = ", sd(pe.vals[[1]] - pe.vals[[2]]), "\nSD of lc.ge.med - lc.ce.med = ", sd(pe.vals[[3]] - pe.vals[[4]]), "\n" )
	}

comb.samps <- function(samp1, samp2)
	{
	uc.ge.med <- c(samp1[[1]], samp2[[1]])
	uc.ce.med <- c(samp1[[2]], samp2[[2]])
	lc.ge.med <- c(samp1[[3]], samp2[[3]])
	lc.ce.med <- c(samp1[[4]], samp2[[4]])
	return(list(uc.ge.med, uc.ce.med, lc.ge.med, lc.ce.med) ) 
	}
	
tum.samps <- read.table("500.tumorsamps")

plot.hists <- function(data = tum.samps)
	{
	uc.n.48 <- data[[1]]
	uc.n.96 <- data[[2]]
	lc.n.48 <- data[[3]]
	lc.n.96 <- data[[4]]
	par(mfrow = c(2,2) )
	hist(lc.n.48, main = c("LC NPMLE", "n = 48"), xlab = "Sampled Estimate" )
	hist(lc.n.96, main = c("LC NPMLE", "n = 96"), xlab = "Sampled Estimate" )
	hist(uc.n.48, main = c("Unconstrained NPMLE", "n = 48"), xlab = "Sampled Estimate" )
	hist(uc.n.96, main = c("Unconstrained", "n = 96"), xlab = "Sampled Estimate" )

	}