#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>

using namespace std;

inline double soft_threshold(double z, double lambda)
{
    if (z > lambda)
    {
        z = z - lambda;
    }
    else if (z < -lambda)
    {
        z = z + lambda;
    }
    else
    {
        z = 0;
    }
    return z;
}

//' coordinate descent algorithm for lasso
//'
//' @param X predictor matrix
//' @param y response vector
//' @param lambda penalty parameter
//' @param max_iteration maximum iteration number
//' @param tol tolanrance parameter
//' @return hat_beta lasso estimation with lambda
//' @export
//' @aliases codelasso
// [[Rcpp::export]]
arma::vec codelasso(arma::mat &X, arma::vec &y, double lambda, unsigned max_iteration, double tol)
{
    arma::uword p = X.n_cols;
    arma::uword n = X.n_rows;
    // initilize hat_beta
    arma::vec hat_beta = arma::zeros(p);
    // initilize residuals
    arma::vec r = y;
    // calculate the norm2 of each column of X (scaled by n) and X'y/n 
    arma::vec X_col_norm2 = arma::zeros(p); 
    for (unsigned j = 0; j < p; j++)
    {
        X_col_norm2(j) = dot(X.col(j), X.col(j)) / n; 
    } 

    unsigned t = 0;
    double gap = INFINITY;
    arma::vec hat_beta_old = arma::zeros(p);
    while ((t < max_iteration) & (gap > tol))
    {
        hat_beta_old = hat_beta;
        // update hat_beta
        for (unsigned j = 0; j < p; j++)
        {
            double temp = dot(X.col(j), r) / n + X_col_norm2(j) * hat_beta(j);
            hat_beta(j) = soft_threshold(temp, lambda) / X_col_norm2(j); 
            // update residuals
            r = y - X * hat_beta;
        }
        // update gap and index t
        gap = norm((hat_beta - hat_beta_old)) / norm(hat_beta_old);
        t = t + 1;
    }
    return hat_beta;
}

//' lasso solution path via coordinate descent algorithm
//' @param X predictor matrix
//' @param y response vector
//' @param nlambda the number of penalty parameter
//' @param max_iteration maximum iteration number
//' @param tol tolanrance parameter
//' @param min_lambda minimum lambda
//' @return A list consists of
//' \item{beta_path}{the solution path of beta}
//' \item{lambda_path}{the lambda path}
//' \item{sigma2_path}{the path of estimated variance via maximum likelihood estimation}
//' @export
//' @aliases codelasso_path
// [[Rcpp::export]]
Rcpp::List codelasso_path(arma::mat &X, arma::vec &y, unsigned nlambda, unsigned max_iteration, double tol, double min_lambda)
{
    unsigned n = X.n_rows;
    unsigned p = X.n_cols;
    double max_lambda = arma::max(abs(X.t() * y)) / n;
    arma::vec Lambda = arma::zeros(nlambda);
    if (min_lambda == 0)
    {
        Lambda.subvec(1, nlambda - 1) = arma::linspace(log(max_lambda / nlambda), log(max_lambda), nlambda - 1);
        Lambda.subvec(1, nlambda - 1) = arma::exp(Lambda.subvec(1, nlambda - 1));
    }
    else
    {
        Lambda = arma::linspace(log(min_lambda), log(max_lambda), nlambda);
        Lambda = arma::exp(Lambda);
    }
    arma::mat beta_path = arma::zeros(p, nlambda);
    for (unsigned k = 0; k < nlambda; k++)
    {
        beta_path.col(k) = codelasso(X, y, Lambda(k), max_iteration, tol);
    }

    arma::vec sigma2_path = arma::zeros(nlambda);
    for (arma::uword k = 0; k < nlambda; k++)
    {
        arma::vec res = y - X * beta_path.col(k);
        sigma2_path(k) = arma::dot(res, res) / n;
    }

    Rcpp::List result;
    result["beta_path"] = beta_path;
    result["lambda_path"] = Lambda;
    result["sigma2_path"] = sigma2_path;
    return result;
}