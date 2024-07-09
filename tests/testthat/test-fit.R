test_that("basic input checks work", {
  set.seed(10)
  n <- 100
  X <- runif(n)
  Y <- (X - 0.5)^2 + 0.5 + rnorm(n, sd = 1e-1)
  model <- "linear"
  family <- "beta"
  
  expect_error(glam("A", Y, model, family))
  expect_error(glam(X, Y, model, "poisson"))
  expect_error(glam(X, Y, model, "normal"))
  expect_error(glam(X, Y, "additive", "binomial"))
})

test_that("glm's capture true coefficients", {
  set.seed(10)
  n <- 2000
  d <- 2
  X <- matrix(runif(n*d), ncol = d)
  beta_true <- c(-1, 2, -1)
  eta_true <- cbind(rep(1, n), X) %*% beta_true
  
  # beta
  mu_true <- 1 / (1 + exp(-eta_true))
  expect_equal(glam(X, mu_true, "linear", "beta")$coef, beta_true)
  
  # binomial
  Y <- rbinom(n, 1, mu_true)
  expect_equal(glam(X, Y, "linear", "binomial")$coef, beta_true, tolerance = 1e-1)
  
  # poisson
  Y <- rpois(n, exp(eta_true))
  expect_equal(glam(X, Y, "linear", "poisson")$coef, beta_true, tolerance = 1e-1)
  
  # gamma
  Y <- exp(eta_true)
  expect_equal(glam(X, Y, "linear", "gamma")$coef, beta_true)
})

test_that("gam's fit as expected", {
  set.seed(10)
  n <- 2000
  d <- 2
  X <- matrix(runif(n*d), ncol = d)
  beta_true <- c(-1, 2, -1)
  eta_true <- c(cbind(rep(1, n), X) %*% beta_true)
  
  # beta
  mu_true <- c(1 / (1 + exp(-eta_true)))
  Y <- mu_true
  expect_equal(glam(X, Y, "additive", "beta")$mu, mu_true, tolerance = 1e-1)
  
  # poisson
  mu_true <- exp(eta_true)
  Y <- rpois(n, mu_true)
  expect_equal(glam(X, Y, "additive", "poisson")$mu, mu_true, tolerance = 1e-1)
  
  # gamma
  mu_true <- exp(eta_true)
  Y <- mu_true
  expect_equal(glam(X, Y, "additive", "gamma")$mu, mu_true, tolerance = 1e-1)
})




