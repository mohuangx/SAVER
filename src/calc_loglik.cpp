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

/*
 calc.loglik.a <- function(a, y, mu, sf) {
  n <- length(y)
  t1 <- 1 / a
  func1 <- n*t1*log(t1)
  func3 <- -n*gammaln(t1)
  
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  
  func2 <- -t1*sum(log(mu))
  func4 <- sum(gammaln(y+t1))
  
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func5 <- -sum((y+t1)*log(sf+t1/mu))
  
  return(-sum(func1, func2, func3, func4, func5))
}
 */

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
double calc_loglik_b(double b, NumericVector y, NumericVector mu, NumericVector sf) {
    double t1 = 1.0 / b;
    double log_t1 = std::log(t1);
    int n = y.length();
    double func1 = 0;
    double func2 = 0;
    double func3 = 0;
    double func4 = 0;
    if (mu.length() == 1) {
        double t2 = mu[0] * t1;
        func1 = t2 * log_t1 * n;
        func2 -= lgamma(t2) * n;
        double log_t_sf = log(sf[0]+t1);
        for (int i = 0; i < n; i++) {
            double t3 = y[i] + t2;
            func3 += lgamma(t3);
            if (sf.length() > 1) {
                func4 -= t3*std::log(sf[i]+t1);
            } else {
                func4 -= t3*log_t_sf;
            }
        }
    } else {
        double log_t_sf = log(sf[0]+t1);
        for (int i = 0; i < n; i++) {
            double t2 = mu[i] * t1;
            func1 += t2 * log_t1;
            func2 -= lgamma(t2);
            double t3 = y[i] + t2;
            func3 += lgamma(t3);
            if (sf.length() > 1) {
                func4 -= t3*std::log(sf[i]+t1);
            } else {
                func4 -= t3*log_t_sf;
            }
        }
    }
    return -(func1+func2+func3+func4);
}

/*
 calc.loglik.k <- function(k, y, mu, sf) {
  t1 <- 1 / k
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  
  func3 <- sum(mu^2*log(mu)*t1)
  func4 <- sum(mu^2*log(t1)*t1)
  func5 <- -sum(lgamma(mu^2*t1))
  
  func6 <- sum(lgamma(y+mu^2*t1))
  
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func7 <- -sum((y+mu^2*t1)*log(sf+mu*t1))
  
  return(-sum(func3, func4, func5, func6, func7))
}
 */

// [[Rcpp::export]]
double calc_loglik_k(double k, NumericVector y, NumericVector mu, NumericVector sf) {
    int n = y.length();
    double t1 = 1.0 / k;
    double log_t1 = std::log(t1);
    
    double func3 = 0;
    double func4 = 0;
    double func5 = 0;
    double func6 = 0;
    double func7 = 0;
    
    if (mu.length() == 1) {
        double mu_square = mu[0] * mu[0];
        double t3 = mu_square * t1;
        double t2 = t3 * n;
        func3 = t2 * log(mu[0]);
        func4 = t2 * log_t1;
        func5 = -lgamma(t3) * n;
        double t4 = log(sf[0]+t1 * mu[0]);
        for (int i = 0; i < n; i++) {
            double y_plus_t3 = y[i] + t3;
            func6 += lgamma(y_plus_t3);
            if (sf.length() > 1) {
                func7 -= y_plus_t3 * log(sf[i]+t1 * mu[0]);
            } else {
                func7 -= y_plus_t3 * t4;
            }
        }
    } else {
        for (int i = 0; i < n; i++) {
            double mu_square = mu[i] * mu[i];
            double t3 = mu_square * t1;
            func3 += t3 * log(mu[i]);
            func4 += t3 * log_t1;
            func5 -= lgamma(t3);
            double y_plus_t3 = y[i] + t3;
            func6 += lgamma(y_plus_t3);
            if (sf.length() > 1) {
                func7 -= y_plus_t3 * log(sf[i]+t1 * mu[i]);
            } else {
                func7 -= y_plus_t3 * log(sf[0]+t1 * mu[i]);
            }
        }
    }
    return -(func3+func4+func5+func6+func7);
}