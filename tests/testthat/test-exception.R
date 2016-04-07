context("Low WAFOM Point Set mvnorm test exception")
library(lwafomps)

test_that("mvnorm stops when arguments are incorrect", {
  param.s <- 4
  param.lower <- rep(-Inf, param.s)
  param.upper <- rep(0.0, param.s)
  param.mean <- rep(0.0, param.s - 1)
  param.shift <- 0
  eps = 0.2 / param.s
  coval<- matrix(eps, nrow=param.s, ncol=param.s)
  d<-diag(1.0 - eps, param.s)
  param.coval <- coval + d
  expect_error(mvnorm(param.s, param.lower, param.upper, param.mean, param.coval),
  "dimension mismatch!", fixed=TRUE)

  param.mean <- rep(0.0, param.s)
  eps = 0.5
  coval<- matrix(eps, nrow=param.s, ncol=param.s)
  d<-diag(1.0 - eps, param.s)
  param.coval <- coval + d
  expect_error(mvnorm(param.s, param.lower, param.upper, param.mean, param.coval),
  "coval is not diagonally dominant", fixed=TRUE)

  param.mean <- rep(0.0, param.s)
  eps = 0.2 / param.s
  coval<- matrix(eps, nrow=param.s, ncol=param.s)
  coval[1,2] <- 0.003
  d<-diag(1.0 - eps, param.s)
  param.coval <- coval + d
  expect_error(mvnorm(param.s, param.lower, param.upper, param.mean, param.coval),
  "coval is not symmetric", fixed=TRUE)

})
