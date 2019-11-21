#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>

using namespace std;

int sgn(double x)
{
    if (x > 0)
    {
        return 1;
    }
    else if (x == 0)
    {
        return 0;
    }
    else
    {
        return -1;
    }
}

//' stagewise lasso 
//' 
//' @param X predictor matrix 
//' @param y response vector 
//' @param eps stepsize 
//' @param tol tolanrance paramater 
//' @param max_iteration maximum iteration number 
//' @return A list consists of 
//' \item{beta_path}{the solution path of beta}
//' \item{lambda_path}{the path of lambda} 
//' \item{sigma2_path}{the path of estimated variance via maximum likelihood estimation}
//' @export 
//' @aliases blasso 
// [[Rcpp::export]]
Rcpp::List blasso(arma::mat &X, arma::mat &y, double eps, double tol, unsigned max_iteration)
{
    unsigned n = X.n_rows;
    unsigned p = X.n_cols;
    arma::mat Beta = arma::zeros(p, max_iteration + 1);
    arma::vec Lambda = arma::zeros(max_iteration + 1);
    arma::uvec active_set = arma::zeros<arma::uvec>(p + 1);
    unsigned active_set_size;

    // initialize residuals
    arma::vec e = y;
    // calculate correlation
    arma::vec corr = X.t() * e;
    // search max correlation
    arma::uword hat_j = index_max(abs(corr));
    // initialize active set
    active_set(0) = hat_j;
    active_set_size = 1;
    // initialize beta0
    double s = sgn(corr(hat_j)) * eps; // set stepsize
    Beta(hat_j, 0) = s;
    // update residuals
    arma::uvec t = arma::zeros<arma::uvec>(1);
    e = e - X.col(hat_j) * Beta(hat_j, 0);
    // initialize lambda
    double delta_loss = s * dot(X.col(hat_j), e) / n - pow(eps, 2) * dot(X.col(hat_j), X.col(hat_j)) / (2 * n);
    Lambda(t(0)) = delta_loss / eps;

    while ((delta_loss > tol) & (t(0) < max_iteration))
    {
        // try backward step
        delta_loss = -INFINITY;
        double delta_loss_temp; 
        double s_hat = 0;
        unsigned j; 
        unsigned hat_j_index = 0; 
        // cout << "t: " << t(0) << endl;
        // search max delta_loss and corresponding index j in active_set
        for (unsigned k = 0; k < active_set_size; k++)
        {
            j = active_set(k);
            s = -sgn(Beta(j, t(0))) * eps;
            delta_loss_temp = s * dot(X.col(j), e) / n - pow(eps, 2) * dot(X.col(j), X.col(j)) / (2 * n);
            if (delta_loss_temp > (delta_loss + tol))
            {
                // update delta_loss and hat_j
                delta_loss = delta_loss_temp;
                hat_j = j;
                s_hat = s;
                hat_j_index = k;
            }
        }
        if (delta_loss > tol)
        {
            // update beta and e
            Beta.col(t(0) + 1) = Beta.col(t(0));
            Beta(hat_j, t(0) + 1) = Beta(hat_j, t(0)) + s_hat;
            e = e - X.col(hat_j) * s_hat;
            // delete hat_j in the active_set
            if (abs(Beta(hat_j, t(0))) < tol)
            {
                active_set.subvec(hat_j_index, active_set_size - 1) = active_set.subvec(hat_j_index + 1, active_set_size);
                active_set_size = active_set_size - 1;
            }
            // update lambda
            Lambda(t(0) + 1) = Lambda(t(0));
            t(0) = t(0) + 1;
        }
        else // forward step
        {
            // search max delta_loss and corresponding index in all set j=1,...,p
            for (unsigned j = 0; j < p; j++)
            {
                delta_loss_temp = eps * abs(dot(X.col(j), e)) / n - pow(eps, 2) * dot(X.col(j), X.col(j)) / (2 * n);
                if (delta_loss_temp > (delta_loss + tol))
                {
                    // update delta_loss and hat_j
                    delta_loss = delta_loss_temp;
                    hat_j = j;
                    s_hat = eps * sgn(dot(X.col(j), e));
                }
            }
            if (delta_loss > tol)
            {
                // update beta and e
                Beta.col(t(0) + 1) = Beta.col(t(0));
                Beta(hat_j, t(0) + 1) = Beta(hat_j, t(0)) + s_hat;
                e = e - X.col(hat_j) * s_hat;
                // update active set
                bool add_j = true;
                for (unsigned k = 0; k < active_set_size; k++)
                {
                    if (active_set(k) == hat_j)
                    {
                        add_j = false;
                        break;
                    }
                }
                if (add_j == true)
                {
                    active_set(active_set_size) = hat_j;
                    active_set_size = active_set_size + 1;
                }
                // update lambda
                Lambda(t(0) + 1) = ((delta_loss - tol) / eps) < Lambda(t(0)) ? ((delta_loss - tol) / eps) : Lambda(t(0));
                t(0) = t(0) + 1;
            }
        }
    }
    arma::mat beta_path = Beta.cols(0, t(0));
    arma::vec lambda_path = Lambda.subvec(0, t(0)); 
    arma::vec sigma2_path = arma::zeros(t(0) + 1); 
    for (arma::uword k = 0; k < t(0) + 1; k++)
    {
        arma::vec res = y - X * beta_path.col(k); 
        sigma2_path(k) = arma::dot(res, res) / n;
    }
    

    Rcpp::List result; 
    result["beta_path"] = beta_path; 
    result["lambda_path"] = lambda_path; 
    result["sigma2_path"] = sigma2_path; 

    return result; 
}