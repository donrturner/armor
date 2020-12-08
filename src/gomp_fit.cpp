// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include "math.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Calculate the partial derivative of f with respect to eta
arma::colvec df_deta(const arma::colvec& x, double eta, double b){
    return(b * exp(eta * (-exp(b * x)) + b * x + eta) % (eta * (exp(b * x) - 1) - 1));
}

// Calculate the partial derivative of f with respect to b
arma::colvec df_db(const arma::colvec& x, double eta, double b){
    return(-eta * (-exp(eta * (-exp(b * x)) + b * x + eta)) % (b * x % (eta * exp(b * x) - 1) - 1));
}

// Gompertz density function
arma::colvec gompfun(const arma::colvec& x, const arma::colvec& par){
    return(par(1) * par(0) * exp(par(0) + par(1) * x - par(0) * exp(par(1) * x)));
}

// Gompertz distribution fit using nonlinear least squares (Levenberg-Marquardt algorithm)
// x - a column vector of ages of length n
// y - a column vector of length n whose elements are the number of deaths for a given age
// par - eta and b, the shape and scale parameters for the Gompertz distribution
// eps - stopping criterion, default value 0.0001
// lambda - damping factor, default value 1
// lamup - what lambda is multiplied with given an unsuccessful step
// lamdown - what lambda is divided by given a successful step
// [[Rcpp::export]]
Rcpp::List gomp_fit_c(const arma::colvec& x, const arma::colvec& y, const arma::colvec& parstart,
                      double eps = 0.0001, double lambda = 1, double lamup = 1.1, double lamdown = 1.5, int maxIter = 100){
    arma::colvec par = parstart;
    arma::colvec obj(2, arma::fill::zeros);
    arma::colvec res = y - gompfun(x, par);
    arma::colvec lamvec(2);
    obj(1) = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
    obj(0) = 2 * obj(1);

    // Update is determined by the equation (J^TJ - lambda * I) %*% d = J^T %*% res,
    // where J is our Jacobian matrix. C is our cost factor. A new step is accepted
    // if the corresponding value of C is lower than the value of C calculated using
    // the old parameters

    arma::mat J, JTJ, lammat, g;
    arma::colvec dC, newpar, newres;
    double C, Cnew;

    int i = 0;

    while(obj(0) - obj(1) > eps){
        if(i == maxIter){
            return(".");
        }
        J = arma::join_rows(df_deta(x, par(0), par(1)), df_db(x, par(0), par(1)));
        JTJ = J.t() * J;
        lamvec.fill(lambda);
        lammat = arma::diagmat(lamvec);
        g = JTJ + lammat;
        dC = J.t() * res;
        C = 0.5 * std::inner_product(res.begin(), res.end(), res.begin(), 0.0);

        // Evaluate new parameters and cost
        newpar = par - inv(g) * dC;
        newres = y - gompfun(x, newpar);
        Cnew = 0.5 * std::inner_product(newres.begin(), newres.end(), newres.begin(), 0.0);

        // If our new C is lower, we accept the update and increase step size.
        if(Cnew < C){
            par = newpar;
            res = newres;
            obj(0) = obj(1);
            obj(1) = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
            lambda /= lamdown;
        // If our new C is not lower, we do not accept the update and decrease the step size
        }else{
            lambda *= lamup;
        }
        i++;
    }
    return Rcpp::List::create(Rcpp::Named("par") = par,
                              Rcpp::Named("MSE") = obj(1));
}
