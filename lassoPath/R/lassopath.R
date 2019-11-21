#' generate lasso path 
#' 
#' @param X predictor matrix 
#' @param y response vector 
#' @param control control list consists of all parameters for stagewise and coordinate descent algorithm 
#' @export 
#' @import Rcpp
#' @aliases lasso.path 
#' @useDynLib lassoPath 
lasso.path <- function(X, y, control = path.control()) {
  tol <- control$tol
  max_iteration <- control$max_iteration
  if (control$algo == "blasso") {
    eps <- control$eps
    fit <- blasso(X, y, eps, tol, max_iteration)
  } else if (control$algo == "code") {
    nlambda <- control$nlambda
    min_lambda <- control$min_lambda
    fit <- codelasso_path(X, y, nlambda, max_iteration, tol, min_lambda)
  } else {
    stop("No this kind of algorithm in the package !")
  }
  Beta <- fit$beta_path
  Lambda <- fit$lambda_path

  return(list(Beta = Beta, Lambda = Lambda))
}

#' configure control list for lasso path 
#' 
#' @param algo algorithm type where "blasso" and "code" mean stagewise and coordinate descent algorithm, respectively.
#' @param eps stepsize for stagewise algorithm 
#' @param tol tolanrance parameter 
#' @param max_iteration maximum iteration number 
#' @param nlambda the number of penalty parameter for coordinate descent algorithm 
#' @param min_lambda minimum penalty parameter for coordinate descent algorithm 
#' @return A list consists of 
#' \item{algo}{algorithm type}
#' \item{eps}{stepsize for stagewise algorithm}
#' \item{tol}{tolanrance parameter} 
#' \item{max_iteration}{maximum iteration number}
#' \item{nlambda}{the number of penalty parameter for coordinate descent algorithm}
#' \item{min_lambda}{minimum penalty parameter for coordinate descent algorithm}
#' @export 
#' @aliases path.control 
path.control <- function(
    algo = "blasso",
    eps = 1e-2,
    tol = 1e-6,
    max_iteration = 5000,
    nlambda = 100,
    min_lambda = 0
) {
  ret <- list(
        algo = algo,
        eps = eps,
        tol = tol,
        max_iteration = max_iteration,
        nlambda = nlambda,
        min_lambda = min_lambda
    )

  return(ret)
} 

.onUnload <- function (libpath) {library.dynam.unload("lassoPath", libpath)}