/*
 * Copyright (C) 2019 Quanli Wang
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 */ 

/*
#' Calculates marginal likelihood
#'
#' Calculates the marginal likelihood given the prediction under constant
#' coefficient of variation (a), Fano factor (b), and variance (k).
#'
#' \code{calc.loglik.a} returns the shifted negative log-likelihood under
#' constant coefficient of variation.
#' \code{calc.loglik.b} returns the shifted negative log-likelihood under
#' constant Fano factor.
#' \code{calc.loglik.k} returns the shifted negative log-likelihood under
#' constant variance.
#'
#' @param a,b,k Prior parameter.
#'
#' @param y A vector of observed gene counts.
#'
#' @param mu A vector of predictions from \code{\link{expr.predict}}.
#'
#' @param sf Vector of normalized size factors.
#'
#' @return A shifted negative marginal log-likelihood.
#'
*/

#include "lgamma.h"

// [[Rcpp::export]]
double calc_loglik_a(double a, NumericVector y, NumericVector mu, NumericVector sf) {
    int n = y.length();
    double t1 = 1.0 / a;
    double log_t1 = std::log(t1);
    
    double func1 = 0;
    double func2 = 0;
    double func3 = 0;
    double func4 = 0;
    double func5 = 0;
    
    func1 = t1 * log_t1 * n;
    func3 = -lgamma(t1) * n;
    
    if (mu.length() == 1) {
        func2 = -log(mu[0]) * t1 * n;
        double t2 = log(sf[0]+t1 / mu[0]);
        for (int i = 0; i < n; i++) {
            double y_plus_t1 = y[i] + t1;
            func4 += lgamma(y_plus_t1);
            if (sf.length() > 1) {
                func5 -= y_plus_t1 * log(sf[i]+t1 / mu[0]);
            } else {
                func5 -= y_plus_t1 * t2;
            }
        }
    } else {
        for (int i = 0; i < n; i++) {
            func2 -= log(mu[i]);
            double y_plus_t1 = y[i] + t1;
            func4 += lgamma(y_plus_t1);
            if (sf.length() > 1) {
                func5 -= y_plus_t1 * log(sf[i]+t1 / mu[i]);
            } else {
                func5 -= y_plus_t1 * log(sf[0]+t1 / mu[i]);
            }
        }
        func2 *= t1;
        
    }
    return -(func1+func2+func3+func4+func5);
}



// [[Rcpp::export]]
NumericVector calc_loglik_b(NumericVector b, NumericVector y, NumericVector mu, NumericVector sf) {
    int n = y.length();
    int n_b = b.length();
    NumericVector log_t1(n_b);
    NumericVector t1(n_b);
    
    for (int i = 0; i < n_b; i++) {
        t1[i] = 1.0 / b[i];
        log_t1[i] =  std::log(t1[i]);
    }
    NumericVector result(n_b);
    
    if (mu.length() == 1) {
        for (int b_index = 0; b_index < n_b; b_index++) {
            double func1 = 0;
            double func2 = 0;
            double func3 = 0;
            double func4 = 0;
            double t2 = mu[0] * t1[b_index];
            func1 = t2 * log_t1[b_index] * n;
            func2 -= lgamma(t2) * n;
            double log_t_sf = log(sf[0]+t1[b_index]);
            for (int i = 0; i < n; i++) {
                double t3 = y[i] + t2;
                func3 += lgamma(t3);
                if (sf.length() > 1) {
                    func4 -= t3*std::log(sf[i]+t1[b_index]);
                } else {
                    func4 -= t3*log_t_sf;
                }
            }
            result[b_index] = -(func1+func2+func3+func4);
        }
    } else {
        for (int b_index = 0; b_index < n_b; b_index++) {
            double log_t_sf = log(sf[0]+t1[b_index]);
            double func1 = 0;
            double func2 = 0;
            double func3 = 0;
            double func4 = 0;
            for (int i = 0; i < n; i++) {
                double t2 = mu[i] * t1[b_index];
                func1 += t2 * log_t1[b_index];
                func2 -= lgamma(t2);
                double t3 = y[i] + t2;
                func3 += lgamma(t3);
                if (sf.length() > 1) {
                    func4 -= t3*std::log(sf[i]+t1[b_index]);
                } else {
                    func4 -= t3*log_t_sf;
                }
            }
            result[b_index] = -(func1+func2+func3+func4);
        }
    }
    return result;
}

// [[Rcpp::export]]
NumericVector calc_loglik_k(NumericVector k, NumericVector y, NumericVector mu, NumericVector sf) {
    int n = y.length();
    
    int n_k = k.length();
    NumericVector t1(n_k);
    NumericVector log_t1(n_k);
    for (int i = 0; i < n_k; i++) {
        t1[i] = 1.0 / k[i];
        log_t1[i] =  std::log(t1[i]);
    }
    NumericVector result(n_k);
    
    if (mu.length() == 1) {
        double mu_square = mu[0] * mu[0];
        double log_mu = log(mu[0]);
        for (int k_index = 0; k_index < n_k; k_index++) {
            double func3 = 0;
            double func4 = 0;
            double func5 = 0;
            double func6 = 0;
            double func7 = 0;
            double t3 = mu_square * t1[k_index];
            double t2 = t3 * n;
            func3 = t2 * log_mu;
            func4 = t2 * log_t1[k_index];
            func5 = -lgamma(t3) * n;
            double t1_mu = t1[k_index] * mu[0];
            double t4 = log(sf[0]+t1_mu);
            for (int i = 0; i < n; i++) {
                double y_plus_t3 = y[i] + t3;
                func6 += lgamma(y_plus_t3);
                if (sf.length() > 1) {
                    func7 -= y_plus_t3 * log(sf[i]+t1_mu);
                } else {
                    func7 -= y_plus_t3 * t4;
                }
            }
            result[k_index] = -(func3+func4+func5+func6+func7);
        }
    } else {
        double* mu_square = new double[n];
        double* log_mu = new double[n];
        for (int i = 0; i < n; i++) {
            mu_square[i] = mu[i] * mu[i];
            log_mu[i] = log(mu[i]);
        }
        for (int k_index = 0; k_index < n_k; k_index++) {
            double func3 = 0;
            double func4 = 0;
            double func5 = 0;
            double func6 = 0;
            double func7 = 0;
            for (int i = 0; i < n; i++) {
                double t3 = mu_square[i] * t1[k_index];
                func3 += t3 * log_mu[i];
                func4 += t3 * log_t1[k_index];
                func5 -= lgamma(t3);
                double y_plus_t3 = y[i] + t3;
                func6 += lgamma(y_plus_t3);
                if (sf.length() > 1) {
                    func7 -= y_plus_t3 * log(sf[i]+t1[k_index] * mu[i]);
                } else {
                    func7 -= y_plus_t3 * log(sf[0]+t1[k_index] * mu[i]);
                }
            }
            result[k_index] = -(func3+func4+func5+func6+func7);
        }
        delete [] mu_square;
        delete [] log_mu;
    }
    return result;
}