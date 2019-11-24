#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

static double threshold = 0.001;

double modify(double& x) {
    return x < threshold ? 0: x;
}

// [[Rcpp::export]]
arma::sp_mat set_zero_if_below(arma::sp_mat x, NumericVector t) {
    threshold = t[0];
    x.transform(modify);
    return x;
}

