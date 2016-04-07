context("Low WAFOM Point Set mvnorm test normal")
library(lwafomps)

test_that("mvnorm normal case 1", {
	param.s <- 4
	param.lower <- rep(-Inf, param.s)
	param.upper <- rep(0, param.s)
	param.mean <- rep(0, param.s)
	param.shift <- 0
	eps = 0.2/param.s
	coval <- matrix(eps, nrow = param.s, ncol = param.s)
	d <- diag(1 - eps, param.s)
	param.coval <- coval + d
	options(digits = 20)
	rs <- mvnorm(param.s, param.lower, param.upper, param.mean, param.coval)
	expect_equal(rs$status, "normal")
	expect_equal(rs$value, expected = 2.05995082121662, tolerance = 1e-06, scale = 1)
	expect_lt(rs$abseps, 1e-06)
	expect_lt(rs$releps, 1e-06)
	expect_equal(rs$confidenceLevel, expected=0.95, tolerance = 1e-06, scale = 1)
})

test_that("mvnorm normal case 2", {
	param.s <- 4
	param.lower <- rep(-Inf, param.s)
	param.upper <- rep(0, param.s)
	param.mean <- rep(0, param.s)
	param.shift <- 0
	eps = 0.2/param.s
	coval <- matrix(eps, nrow = param.s, ncol = param.s)
	d <- diag(1 - eps, param.s)
	param.coval <- coval + d
	options(digits = 20)
	rs <- mvnorm(param.s, param.lower, param.upper, param.mean, param.coval, confidenceLevel=0.98)
	expect_equal(rs$status, "normal")
	expect_equal(rs$value, expected = 2.0599507948821, tolerance = 1e-06, scale = 1)
	expect_lt(rs$abseps, 1e-06)
	expect_lt(rs$releps, 1e-06)
	expect_equal(rs$confidenceLevel, expected=0.98, tolerance = 1e-06, scale = 1)
})
