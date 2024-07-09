#' Generalized Linear/Additive Model (GLAM)
#' 
#' @description
#' `glam` fits Generalized Linear Models (GLMs) and Generalized Additive Models (GAMs).
#' 
#' @import gam
#' @import stats
#' @param x a `matrix` with covariate data.
#' @param y a `vector` of responses. 
#' @param model a string stating whether a GLM (`"linear"`) or GAM (`"additive"`) should be fit.
#' @param family the distributional family from which to model response. Options are 
#' `"beta"`, `"binomial"` (GLM only), `"poisson"`, and `"gamma"`.
#' @param intercept a `logical` indicating whether an intercept is desired. Only applicable when
#' `model = "linear"`. Default is `TRUE`.
#' @param tol the tolerance for convergence. Must be positive. Default is `1e-8`.
#' @param max_iter an `integer` specifying the maximum number of allowed iterations. Default is `100`.
#' @details This is a function for training and fitting Generalized Linear Models (GLMs) and
#' Generalized Additive Models (GAMs). It implements these models using Iterative Reweighted Least Squares
#' (IRLS) described in Hastie and Tibshirani 1990 (<doi:10.1214/ss/1177013604>). This function supports Beta
#' regression, Logistic regression, Poisson regression, and Gamma regression (Logistic GAMs are
#' currently not supported).
#' @return The output is a `list` containing:\tabular{ll}{
#'    \code{eta} \tab a `vector` of un-transformed fitted values. \cr
#'    \tab \cr
#'    \code{mu} \tab a `vector` of transformed fitted values. \cr
#'    \tab \cr
#'    \code{num_iter} \tab the number of iterations until convergence or timeout. \cr
#'    \tab \cr
#'    \code{dev} \tab the convergence criteria achieved at the end. \cr
#'    \code{coef} \tab a numeric vector of estimated coefficients 
#'    (only applicable when `model = "linear"`). \cr
#'    }
#' @export
#' 
#' @examples
#' # generate random inputs
#' set.seed(10)
#' n <- 200
#' X <- sort(runif(n))

#' # Beta GLM vs. GAM
#' Y <- (X - 0.5)^2 + 0.5 + rnorm(n, sd = 1e-1)
#' beta_glm <- glam(cbind(X, X^2), Y, model = "linear", family = "beta")
#' beta_gam <- glam(cbind(X, X^2), Y, model = "additive", family = "beta")
#' plot(X, Y, pch = 20, xlab = "X", ylab = "Y", main = "Beta Regression Example")
#' lines(X, beta_glm$mu, col = "red", lwd = 2)
#' lines(X, beta_gam$mu, col = "blue", lwd = 2)
#' legend("bottomright", legend = c("GLM", "GAM"), col = c("red", "blue"), lwd = 2)
glam <- function(x, y, model, family, intercept = TRUE, tol = 1e-8, max_iter = 100){
  
  # check input types
  if(!is.numeric(x)){
    stop("x must be of type 'numeric'.")
  }
  if(!is.numeric(y)){
    stop("y must be of type 'numeric' (binary response data must be coded accordingly).")
  }
  if(!(model %in% c("linear", "additive"))){
    stop("Only supported models are 'linear' for GLM and 'additive' for GAM.")
  }
  if(!(family %in% c("beta", "binomial", "poisson", "gamma"))){
    stop("Currently supported families are Beta, Binomial, Poisson and Gamma.")
  }
  if(!is.numeric(tol)){
    stop("tol must be of type 'numeric'.")
  }
  if(tol <= 0){
    stop("tol must be a positive number.")
  }
  if(!is.numeric(max_iter)){
    stop("max_iter must be of type 'numeric'.")
  }
  if(max_iter <= 1){
    stop("max_iter must be a positive integer.")
  }
  # remove when logistic gam added as option
  if(model == "additive" & family == "binomial"){
    stop("Logstic GAM not supported in current package version.")
  }
  
  # convert to matrix if x is vector
  if(is.null(nrow(x))) x <- as.matrix(x)
  
  # add intercept to design matrix if desired
  if(model == "additive") intercept = F
  if(intercept) x <- cbind(rep(1, nrow(x)), x)
  
  # beta regression
  if(family == "beta"){
    if(any(y < 0) || any(y > 1)){
      stop("Response for Beta regression should be between 0 and 1 inclusive.")
    }
    g = function(mu) log(mu / (1 - mu))
    gi = function(eta) 1 / (1 + exp(-eta))
    zf = function(eta, mu) eta + ((y - mu) / (mu * (1 - mu)))
    wf = function(eta, mu) mu * (1 - mu) * 2
  }
  # logistic regression
  if(family == "binomial"){
    if(length(unique(y)) != 2){
      stop("Response for Binomial regression should contain two unique values.")
    }
    g = function(mu) log(mu / (1 - mu))
    gi <- function(eta) 1 / (1 + exp(-eta))
    zf = function(eta, mu) eta + ((y - mu) / (mu * (1 - mu)))
    wf = function(eta, mu) mu * (1 - mu)
  }
  # poisson regression
  if(family == "poisson"){
    if(any(y < 0) || !(all(y == round(y)))){
      stop("Response for Poisson regression should contain non-negative counts.")
    }
    g = function(mu) log(mu)
    gi = function(eta) exp(eta)
    zf = function(eta, mu) eta + ((y - mu) / mu)
    wf = function(eta, mu) mu
  }
  # gamma regression
  if(family == "gamma"){
    if(any(y <= 0)){
      stop("Response for Gamma regression should contain positive values.")
    }
    g = function(mu) log(mu)
    gi = function(eta) exp(eta)
    zf = function(eta, mu) eta + ((y - mu) / mu)
    wf = function(eta, mu) rep(1, length(eta))
  }
  
  if(model == "linear"){
    return(advr(x, y, 
                gi = gi,
                zf = zf,
                wf = wf, 
                max_iter = max_iter,
                tol = tol))
  }
  if(model == "additive"){
    return(loc_score(x, y,
                     g = g,
                     gi = gi, 
                     zf = zf,
                     wf = wf,
                     max_iter = max_iter,
                     tol = tol))
  }
}