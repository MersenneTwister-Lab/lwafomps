## ------------------------------------------------------------------------
        library(lwafomps)
        param.s <- 4
        param.lower <- rep(-Inf, param.s)
        param.upper <- rep(0.0, param.s)
        param.mean <- rep(0.0, param.s)
        eps = 0.2 / param.s
        coval<- matrix(eps, nrow=param.s, ncol=param.s)
        d<-diag(1.0 - eps, param.s)
        param.coval <- coval + d
        options(digits=20)
        rs <- mvnorm(dimension=param.s,
                     lower=param.lower,
                     upper=param.upper,
                     mean=param.mean,
                     coval=param.coval)
        rs

## ------------------------------------------------------------------------
        library(lwafomps)
        param.s <- 4
        param.lower <- rep(-Inf, param.s)
        param.upper <- rep(0.0, param.s)
        param.mean <- rep(0.0, param.s)
        eps = 0.2 / param.s
        coval<- matrix(eps, nrow=param.s, ncol=param.s)
        d<-diag(1.0 - eps, param.s)
        param.coval <- coval + d
        options(digits=20)
        rs <- mvnorm(dimension=param.s,
                     lower=param.lower,
                     upper=param.upper,
                     mean=param.mean,
                     coval=param.coval,
                     abseps=1.0e-6)
        rs

## ------------------------------------------------------------------------
        library(lwafomps)
        param.s <- 4
        param.lower <- rep(-Inf, param.s)
        param.upper <- rep(0.0, param.s)
        param.mean <- rep(0.0, param.s)
        eps = 0.2 / param.s
        coval<- matrix(eps, nrow=param.s, ncol=param.s)
        d<-diag(1.0 - eps, param.s)
        param.coval <- coval + d
        options(digits=20)
        rs <- mvnorm(dimension=param.s,
                     lower=param.lower,
                     upper=param.upper,
                     mean=param.mean,
                     coval=param.coval,
                     confidenceLevel=0.98)
        rs

## ------------------------------------------------------------------------
        library(lwafomps)
        param.s <- 4
        param.lower <- rep(-Inf, param.s)
        param.upper <- rep(0.0, param.s)
        param.mean <- rep(0.0, param.s)
        eps = 0.2 / param.s
        coval<- matrix(eps, nrow=param.s, ncol=param.s)
        d<-diag(1.0 - eps, param.s)
        param.coval <- coval + d
        options(digits=20)
        rs <- mvnorm(dimension=param.s,
                     lower=param.lower,
                     upper=param.upper,
                     mean=param.mean,
                     coval=param.coval,
                     confidenceLevel=0.98,
                     timeLimit=1.5)
        rs

## ------------------------------------------------------------------------
        library(lwafomps)
        param.s <- 4
        param.lower <- rep(-Inf, param.s)
        param.upper <- rep(0.0, param.s)
        param.mean <- rep(0.0, param.s)
        eps = 0.2 / param.s
        coval<- matrix(eps, nrow=param.s, ncol=param.s)
        d<-diag(1.0 - eps, param.s)
        param.coval <- coval + d
        options(digits=20)
        rs <- mvnorm(dimension=param.s,
                     lower=param.lower,
                     upper=param.upper,
                     mean=param.mean,
                     coval=param.coval,
                     confidenceLevel=0.98,
                     parallelism=32)
        rs

