/*
#' Optimizes variance
#'
#' Finds the prior parameter that maximizes the marginal likelihood given
#' the prediction.
#'
#' \code{calc.a} returns a prior alpha parameter assuming constant
#' coefficient of variation. \code{calc.b} returns a prior beta parameter
#' assuming constant Fano factor. \code{calc.k} returns a prior variance
#' parameter assuming constant variance.
#'
#' @param y A vector of observed gene counts.
#'
#' @param mu A vector of predictions from \code{\link{expr.predict}}.
#'
#' @param sf Vector of normalized size factors.
#'
#' @return A vector with the optimized parameter and the negative
#' log-likelihood.
#'
#'
#' @importFrom  stats optimize ppoints uniroot var
#'
#' @rdname optimize_variance
#' @export
*/

#include "lgamma.h"

typedef NumericVector (*funcPtr)(NumericVector, NumericVector, NumericVector, NumericVector);
typedef NumericVector (*funcShiftPtr)(NumericVector, NumericVector, NumericVector, NumericVector, double);

NumericVector calc_loglik_a_shift(NumericVector a, NumericVector y, NumericVector mu, NumericVector sf, double shift) {
    return calc_loglik_a(a, y, mu, sf) + shift;
}

NumericVector calc_loglik_b_shift(NumericVector a, NumericVector y, NumericVector mu, NumericVector sf, double shift) {
    return calc_loglik_b(a, y, mu, sf) + shift;
}

NumericVector calc_loglik_k_shift(NumericVector a, NumericVector y, NumericVector mu, NumericVector sf, double shift) {
    return calc_loglik_k(a, y, mu, sf) + shift;
}

// [[Rcpp::export]]
NumericVector calc_abk(NumericVector y, NumericVector mu, NumericVector sf, std::string type, NumericVector samp) {
    NumericVector result(2);
    int n = y.length();
    double sum = 0;
    for (int i = 0; i < mu.length();i++) {
        sum += mu[i];
    }
    if (sum == 0 || n ==1) {
        return result;
    }

    double interval_max = 0;
    double * y_over_sf = new double[n];
    double mean = 0;
    for (int i = 0; i < n; i++) {
        if (sf.length() > 0) {
            y_over_sf[i] = y[i] / sf[i];
        } else {
             y_over_sf[i] = y[i] / sf[0];
        }
        mean += y_over_sf[i];
    }
    mean /= n;
    
    double var = 0;
    for (int i = 0; i < n; i++) {
        var += (y_over_sf[i] - mean) * (y_over_sf[i] - mean);
    }
    var /= (n-1);
    delete [] y_over_sf;
    
    funcPtr f = NULL;
    funcShiftPtr f_Shift = NULL;
    if (type == "a") {
      interval_max  = var/mean/mean;
      f = calc_loglik_a;
      f_Shift = calc_loglik_a_shift;
    } else if (type == "b") {
      interval_max = var/mean;
      f = calc_loglik_b;
      f_Shift = calc_loglik_b_shift;
    } else {
      interval_max = var;
      f = calc_loglik_k;
      f_Shift = calc_loglik_k_shift;
    }
    Environment stats("package:stats"); 
    Function optimize = stats["optimize"];
    NumericVector interval(2);
    interval[0] = 0; interval[1] = interval_max;
    List v_vec = optimize(Rcpp::_["f"] = Rcpp::InternalFunction(f), Rcpp::_["interval"] =  interval, y, mu, sf);
    double v = v_vec["minimum"];
    double loglik = v_vec["objective"];

    NumericVector x(2); x[0] = 1e-05; x[1] = v;
    NumericVector r = -f(x,y,mu,sf);
    double min_v = r[0];
    double mle_v = r[1];
    
    if (mle_v - min_v < 0.5) {
        NumericVector xx(1); xx[0] = interval_max;
        NumericVector r_max = -f(xx,y,mu,sf);
        double max_v = r_max[0];
        double v_max = interval_max;
        if (mle_v- max_v > 10.0) {
            Function uniroot = stats["uniroot"];
            interval[0] = 1e-05;
            List uniroot_result = uniroot(Rcpp::_["f"] = Rcpp::InternalFunction(f_Shift), Rcpp::_["interval"] =  interval, y, mu, sf, mle_v - 10);
            v_max = uniroot_result["root"];
        } 
        NumericVector samp_scaled = samp * v_max;
        NumericVector loglik_samp = f(samp_scaled, y, mu, sf);
        double min = -loglik_samp[0];
        for (int i = 1; i < loglik_samp.length();i++) {
            if (-loglik_samp[i] < min) {
                min = -loglik_samp[i];
            }
        }
        double sum_loglik2 = 0;
        v = 0;
        for (int i = 0; i < loglik_samp.length();i++) {
            double loglik2 = exp(-loglik_samp[i] - min);
            v += samp_scaled[i] * loglik2;
            sum_loglik2 += loglik2;
        }
        v /= sum_loglik2;
        NumericVector vv(1); vv[0] = v;
        loglik = f(vv, y, mu, sf)[0];
    }
    if (type == "a" || type == "b") {
        result[0] = 1.0/ v;
    } else {
        result[0] = v;
    }
    result[1] = loglik;
    return result;
}
