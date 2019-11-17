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


#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

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
            if (sf.length() > 0) {
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
            if (sf.length() > 0) {
                func4 -= t3*std::log(sf[i]+t1);
            } else {
                func4 -= t3*log_t_sf;
            }
        }
    }
    return -(func1+func2+func3+func4);
}
