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

##' Multivariate Normal Distribution with Low WAFOM Point Set
##'
##' Compute the distribution function of the multivariate normal
##' distribution for arbitrary limits and correlation matrices.
##'
##'@param lower the vector of lower limits.
##'@param upper the vector of upper limits.
##'@param mean the mean vector.
##'@param coval the corelation matrix.
##'@param abseps expected absolute error tolerance.
##'@param releps expected relative error tolerance.
##'@param probability the expected probability that the error of return value
##'       is smaller than the expected error tolerance.
##'@param tlimit time limit (seconds).
##'@param parallelism degree of parallelism for fast calculation.
##'@return integrated value of multivariate normal distribution.
##'@export
mvnorm <- function(lower, upper, mean, coval, abseps,
                   releps, probability=0.9545, tlimit=Inf, parallelism=4)
{
    start = proc.time()
    count <- 4
    nameList <- c("value", "abseps", "releps", "probability", "status")
    randoms <- sample(2^31, count, replace=TRUE) -1
    pre <- rcppPrecompute(lower, upper, mean, coval)
    result <- foreach(i = 1:count, .combine="c") %do% {
        rcppMvnorm(pre, randoms[i])
    }
    m <- mean(result)
    sdv <- sd(result)
    if (probability > 0.9545) {
        p <- 0.9973
        ae <- 3 * sdv
    } else if (probability > 0.6827) {
        p <- 0.9545
        ae <- 2 * sdv
    } else {
        p <- 0.6827
        ae <- sdv
    }
    re <- ae / m
    if (missing(abseps) && missing(releps)) {
        r <- list(m, ae, re, p, "normal")
        names(r) <- nameList
        return(r)
    }
    if (missing(abseps)) {
        abseps <- m * releps
    }
    if (ae <= abseps ) {
        r <- c(m, ae, re, p, "normal")
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
    while (proc.time() - start < tlimit) {
        m <- mean(result)
        sdv <- sd(result)
        if (probability >= 0.99) {
            p <- 0.9973
            ae <- 3 * sdv
        } else if (probability >= 0.95) {
            p <- 0.9545
            ae <- 2 * sdv
        } else {
            p <- 0.6827
            ae <- sdv
        }
        re <- ae / m
        if (ae <= abseps) {
            r <- list(m, ae, re, p, "normal")
            names(r) <- nameList
            return(r)
        }
        if (m == oldm && sdv == oldsdv) {
            r <- list(m, ae, re, p, "convergence")
            names(r) <- nameList
            return(r)
        }
        total <- total + count
        if (total >= 2^20) {
            r <- list(m, ae, re, p, "total count over")
            names(r) <- nameList
            return(r)
        }
        oldm <- m
        oldsdv <- sdv
        randoms <- sample(2^31, count, replace=TRUE) -1
        result <- append(result, foreach(i = 1:count, .combine="c") %do% {
            rcppMvnorm(pre, randoms[i])
        })
    }
    r <- list(m, ae, re, p, "timeout")
    names(r) <- nameList
    return(r)
}
