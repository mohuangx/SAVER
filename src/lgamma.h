#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

double lgamma(double x);
NumericVector calc_loglik_a(NumericVector a, NumericVector y, NumericVector mu, NumericVector sf);
NumericVector calc_loglik_b(NumericVector a, NumericVector y, NumericVector mu, NumericVector sf);
NumericVector calc_loglik_k(NumericVector a, NumericVector y, NumericVector mu, NumericVector sf);

