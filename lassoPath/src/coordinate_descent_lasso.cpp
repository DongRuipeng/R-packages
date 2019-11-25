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
//' @param hat_beta a warm start for iteration 
//' @param lambda penalty parameter
//' @param max_iteration maximum iteration number
//' @param tol tolanrance parameter
//' @return hat_beta lasso estimation with lambda
//' @export
//' @aliases codelasso
// [[Rcpp::export]]
arma::vec codelasso(arma::mat &X, arma::vec &y, arma::vec &hat_beta, double lambda, unsigned max_iteration, double tol)
{
    arma::uword p = X.n_cols;
    arma::uword n = X.n_rows;
    // initilize hat_beta, gram matrix and Xy (scaled by n)
    arma::mat gram_matrix = arma::zeros(p, p);
    for (unsigned i = 0; i < p; i++)
    {
        for (unsigned j = i; j < p; j++)
        {
            gram_matrix.at(i, j) = arma::dot(X.col(i), X.col(j)) / n;
            gram_matrix.at(j, i) = gram_matrix.at(i, j);
        }
    }
    arma::vec Xy = X.t() * y / n;
    // initialize active set
    arma::uvec active_set = arma::ones<arma::uvec>(p);

    unsigned t = 0;
    double gap = INFINITY;
    arma::vec hat_beta_old;
    arma::uvec ind_j; 
    arma::uvec active_index; 
    while ((t < max_iteration) & (gap > tol))
    {
        hat_beta_old = hat_beta;
        // update hat_beta
        for (unsigned j = 0; j < p; j++)
        {
            ind_j = {j}; 
            active_index = find(active_set == 1); 
            double temp1 = Xy(j) - dot(gram_matrix.submat(active_index, ind_j), hat_beta(active_index));
            double temp2 = temp1 + gram_matrix.at(j, j) * hat_beta(j);
            hat_beta(j) = soft_threshold(temp2, lambda) / gram_matrix.at(j, j);
            if (abs(hat_beta(j)) < tol)
            {
                active_set(j) = 0;
            }
            else
            {
                active_set(j) = 1;
            }
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
    arma::vec hat_beta;
    arma::solve(hat_beta, X, y, arma::solve_opts::fast);
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
        beta_path.col(k) = codelasso(X, y, hat_beta, Lambda(k), max_iteration, tol);
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