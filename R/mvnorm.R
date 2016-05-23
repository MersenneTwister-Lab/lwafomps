##'Multivariate Normal Distribution with Low WAFOM Point Set
##'
##' Compute the distribution function of the multivariate normal
##' distribution for arbitrary limits and correlation matrices.
##'
##'@name lwafomps-package
##'@aliases lwafomps-package lwafomps
##'@docType package
##'@import Rcpp foreach doSNOW
##'@useDynLib lwafomps
NULL

checkParam <- function(dimension, lower, upper, mean, coval, abseps, releps, confidenceLevel,
	timeLimit, parallelism) {
	if (dimension != length(lower) || dimension != length(upper) || dimension != length(mean) ||
		dimension != ncol(coval) || dimension != nrow(coval)) {
		stop("dimension mismatch!")
	}
	e <- try(chol(coval, pivot=TRUE), silent = TRUE)
	if (class(e) == "try-error") {
		stop("coval is not positive semi-definete")
	}
}

calcConfidenceInterval <- function(samples, confidenceLevel) {
	m <- mean(samples)
	s <- sd(samples)
	a <- (1 - confidenceLevel)/2
	n <- length(samples)
	t <- -qt(a, n - 1)
	e <- t * s/sqrt(n)
	c(m, e)
}

##' Multivariate Normal Distribution with Low WAFOM Point Set
##'
##' Compute the distribution function of the multivariate normal
##' distribution for arbitrary limits and correlation matrices.
##'
##'@param dimension dimension or number of variables.
##'@param lower the vector of lower limits.
##'@param upper the vector of upper limits.
##'@param mean the mean vector.
##'@param coval the corelation matrix.
##'@param abseps expected absolute error tolerance.
##'@param releps expected relative error tolerance.
##'@param confidenceLevel confidence level.
##'@param timeLimit time limit (seconds).
##'@param parallelism degree of parallelism for fast calculation.
##'@return integrated value of multivariate normal distribution.
##'@export
mvnorm <- function(dimension, lower, upper, mean, coval, abseps, releps, confidenceLevel = 0.95,
	timeLimit = Inf, parallelism = 4) {
	start = proc.time()
	checkParam(dimension, lower, upper, mean, coval, abseps, releps, confidenceLevel, timeLimit,
		parallelism)
	count <- 4
	nameList <- c("value", "abseps", "releps", "confidenceLevel", "status")
	randoms <- sample(2^31, count, replace = TRUE) - 1
	pre <- rcppPrecompute(lower, upper, mean, coval)
	result <- foreach(i = 1:count, .combine = "c") %do% {
		rcppMvnorm(pre, randoms[i])
	}
	confidence <- calcConfidenceInterval(result, confidenceLevel)
	m <- confidence[1]
	ae <- confidence[2]
	re <- ae/m
	if (missing(abseps) && missing(releps)) {
		r <- list(m, ae, re, confidenceLevel, "normal")
		names(r) <- nameList
		return(r)
	}
	if (missing(abseps)) {
		abseps <- m * releps
	}
	if (ae <= abseps) {
		r <- c(m, ae, re, confidenceLevel, "normal")
		names(r) <- nameList
		return(r)
	}
	#
	# Main Loop
    #
    count <- parallelism
	total <- 0
	oldm <- 0
	oldsdv <- 0
	while (proc.time() - start < timeLimit) {
		confidence <- calcConfidenceInterval(result, confidenceLevel)
		m <- confidence[1]
		ae <- confidence[2]
		re <- ae/m
		if (ae <= abseps) {
			r <- list(m, ae, re, confidenceLevel, "normal")
			names(r) <- nameList
			return(r)
		}
		if (m == oldm && ae == oldae) {
			r <- list(m, ae, re, confidenceLevel, "convergence")
			names(r) <- nameList
			return(r)
		}
		total <- total + count
		if (total >= 2^20) {
			r <- list(m, ae, re, confidenceLevel, "total count over")
			names(r) <- nameList
			return(r)
		}
		oldm <- m
		oldae <- ae
		randoms <- sample(2^31, count, replace = TRUE) - 1
		result <- append(result, foreach(i = 1:count, .combine = "c") %do% {
			rcppMvnorm(pre, randoms[i])
		})
	}
	r <- list(m, ae, re, confidenceLevel, "timeout")
	names(r) <- nameList
	return(r)
}
