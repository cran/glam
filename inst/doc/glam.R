## ----eval = F-----------------------------------------------------------------
#  glam(X, Y, model, family, intercept, tol, max_iter)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ----setup--------------------------------------------------------------------
library(glam)

## -----------------------------------------------------------------------------
# generate random inputs
set.seed(10)
n <- 50
X <- sort(runif(n))

## ----fig.align = "center", fig.height = 5, fig.width = 5----------------------
Y <- (X - 0.5)^2 + 0.5 + rnorm(n, sd = 1e-1)
beta_glm <- glam(cbind(X, X^2), Y, model = "linear", family = "beta")
beta_gam <- glam(cbind(X, X^2), Y, model = "additive", family = "beta")
plot(X, Y, xlab = "X", ylab = "Y", main = "Beta Regression Example")
lines(X, beta_glm$mu, col = "red", lwd = 2)
lines(X, beta_gam$mu, col = "blue", lwd = 2)
legend("bottomright", legend = c("GLM", "GAM"), col = c("red", "blue"), lwd = 2)

## -----------------------------------------------------------------------------
beta_glm$coef

## ----fig.align = "center", fig.height = 5, fig.width = 5----------------------
Y <- rep(0, n)
Y[which(X <= 1/3)] <- rbinom(sum(X <= 1/3), 1, 0.9)
Y[which(X >= 2/3)] <- rbinom(sum(X >= 2/3), 1, 0.9)
log_glm <- glam(cbind(X, X^2), Y, model = "linear", family = "binomial")
plot(X, Y, xlab = "X", ylab = "Y", main = "Binomial Regression Example")
lines(X, log_glm$mu, col = "red", lwd = 2)
legend("bottomright", legend = c("GLM"), col = c("red"), lwd = 2)

## ----fig.align = "center", fig.height = 5, fig.width = 5----------------------
Y <- rpois(n, 3*rnorm(n, mean = 2*X, sd = 1e-4))
pois_glm <- glam(cbind(X, X^2), Y, model = "linear", family = "poisson")
pois_gam <- glam(cbind(X, X^2), Y, model = "additive", family = "poisson")
plot(X, Y, xlab = "X", ylab = "Y", main = "Poisson Regression Example")
lines(X, pois_glm$mu, col = "red", lwd = 2)
lines(X, pois_gam$mu, col = "blue", lwd = 2)
legend("bottomright", legend = c("GLM", "GAM"), col = c("red", "blue"), lwd = 2)

## ----fig.align = "center", fig.height = 5, fig.width = 5----------------------
Y <- rgamma(n, shape = 10, scale = 5*X)
gamma_glm <- glam(cbind(X, X^2), Y, model = "linear", family = "gamma")
gamma_gam <- glam(cbind(X, X^2), Y, model = "additive", family = "gamma")
plot(X, Y, xlab = "X", ylab = "Y", main = "Gamma Regression Example")
lines(X, gamma_glm$mu, col = "red", lwd = 2)
lines(X, gamma_gam$mu, col = "blue", lwd = 2)
legend("bottomright", legend = c("GLM", "GAM"), col = c("red", "blue"), lwd = 2)

