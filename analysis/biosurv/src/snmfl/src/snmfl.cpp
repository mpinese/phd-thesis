// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

extern "C" int nnls_c(double* a, const int* mda, const int* m, const int* n, double* b, double* x, double* rnorm, double* w, double* zz, int* index, int* mode);
arma::vec nnls(arma::mat A, arma::vec b);


void r_message(std::string text);
void r_warning(std::string text);

void r_message(std::string text)
{
    Rcpp::Function msg("message");
    msg(text);
}

void r_warning(std::string text)
{
    Rcpp::Function warn("warning");
    warn(text);
}

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// SNMF/L matrix factorization
//
// Calculates W, H so that A ~= W*H, A >= 0, W >= 0, H >= 0, and
// rk(W) = rk(H) = k.  eta and alpha control fit regularization.
//     
// A:             Matrix to factorize.
// k:             Rank of the factorization.
// eta:           H norm penalty coefficient.  Negative numbers are
//                automatically set to a reasonable value.
// alpha:         W sparsity penalty coefficient.  Negative numbers are
//                automatically set to a reasonable value.
// max_iter:      Maximum number of fit iterations
// conv_delta:    Convergence criterion.  If the improvement in fit
//                residual ||A - WH||_{F}^2 <= conv_delta, consider
//                the fit to have converged.
// conv_interval: Convergence check interval.  Checks for convergence
//                will only be performed every conv_interval iterations.
//
// Returns a list of W, H, norm, and converged.
//
// [[Rcpp::export]]
Rcpp::List snmfl_atomic(const arma::mat& A, arma::uword k, double eta, double alpha, size_t max_iter, double conv_delta, size_t conv_interval)
{
    size_t i;
    arma::uword j;
    double current_residual_norm, last_residual_norm;
    bool converged;

    arma::uword n = A.n_rows;
    arma::uword p = A.n_cols;

    arma::mat W = arma::randu<arma::mat>(n, k) * A.max();
    arma::mat H = arma::zeros<arma::mat>(k, p);
    arma::mat P = arma::zeros<arma::mat>(n + k, k);
    arma::mat Q = arma::zeros<arma::mat>(p + 1, k);

    // Set eta automatically to a reasonable number if < 0
    if (eta < 0)
        eta = A.max();

    // Set alpha automatically to a reasonable number if < 0
    if (alpha < 0)
        alpha = 0.01;

    // Initialize constant sections of P and Q
    P.rows(n, n + k - 1) = arma::eye<arma::mat>(k, k) * sqrt(eta);
    Q.row(p) = arma::ones<arma::rowvec>(k) * sqrt(alpha);

    converged = false;
    for (i = 0; i <= max_iter; i++)
    {
        if (i == 0 || i % conv_interval == 0 || i == max_iter - 1)
        {
            current_residual_norm = arma::accu(arma::square(A - W*H));
            if (i != 0)
            {
                // Convergence check
                if (last_residual_norm - current_residual_norm <= conv_delta)
                {
                    converged = true;
                    break;
                }
            }

            if (i == max_iter - 1)
                break;
        }
        last_residual_norm = current_residual_norm;

        // Calculate the next H
        P.rows(0, n-1) = W;
        for (j = 0; j < p; j++)
            H.col(j) = nnls(P, arma::join_cols(A.col(j), arma::zeros<arma::vec>(k)));

        // Calculate the next W
        Q.rows(0, p-1) = H.t();
        for (j = 0; j < n; j++)
            W.row(j) = nnls(Q, arma::join_cols(A.row(j).t(), arma::zeros<arma::vec>(1))).t();
    }

    return Rcpp::List::create(
        Rcpp::Named("W") = W,
        Rcpp::Named("H") = H,
        Rcpp::Named("norm") = current_residual_norm,
        Rcpp::Named("converged") = converged);
}


arma::vec nnls(arma::mat A, arma::vec b)
{
    int mda, m, n, mode, i;
    int *index;
    double rnorm;
    double *x;
    double *w;
    double *zz;

    mda = m = A.n_rows;
    n = A.n_cols;

    arma::mat Acopy(A);
    arma::vec bcopy(b);

    arma::vec result(n);
    result.fill(arma::datum::nan);

    rnorm = 0;
    mode = 0;

    if ((x = new double[n]) == NULL)
    {
        Rcpp::stop("nnls: Memory allocation failed");
    }

    if ((w = new double[n]) == NULL)
    {
        delete x;
        Rcpp::stop("nnls: Memory allocation failed");
    }

    if ((index = new int[n]) == NULL)
    {
        delete w;
        delete x;
        Rcpp::stop("nnls: Memory allocation failed");
    }

    if ((zz = new double[m]) == NULL)
    {
        delete index;
        delete w;
        delete x;
        Rcpp::stop("nnls: Memory allocation failed");
    }

    for (i = 0; i < n; i++)
    {
        w[i] = 0;
        x[i] = 0;
        index[i] = 0;
    }

    for (i = 0; i < m; i++)
        zz[i] = 0;

    mda = m = A.n_rows;
    n = A.n_cols;

    nnls_c(Acopy.memptr(), &mda, &m, &n, bcopy.memptr(), x, &rnorm, w, zz, index, &mode);

    for (i = 0; i < n; i++)
        result(i) = x[i];

    delete zz;
    delete index;
    delete w;
    delete x;

    return result;
}
