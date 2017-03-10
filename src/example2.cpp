#include <Rcpp.h>
using namespace Rcpp;

//' Rcpp Example
//'
//' This is an example
//'
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


//' Rcpp Example
//'
//' This is an example
//'
//' @export
// [[Rcpp::export]]
NumericVector dnorm_cpp(NumericVector x, double mean=0, double sd=1) {
    double pi = 3.141592653589793238463;
    return -(log(2 * pi) + 2 * log(sd) + pow((x-mean)/sd,2))/2;
}

//' Array
//'
//' Array in rcpp
//'
//' @export
// [[Rcpp::export]]
NumericVector arrayC(NumericVector input, IntegerVector dim) {
    input.attr("dim") = dim;
    return input;
}


//' Prediction negLogLik Function
//'
//' Prediction negLogLik Function used in in pred function in pred.R
//'
//' @export
// [[Rcpp::export]]
double pnll_cpp(double        fv,
                NumericVector score_row,
                NumericVector wts_row,
                NumericVector err,
                NumericVector mu,
                NumericVector fl) {
    //double pi = 3.141592653589793238463;
    NumericVector nll1 = wts_row*(0.79817986835+2*log(err)+pow((score_row-mu - fl * fv)/err,2))/2;
    // int n = nll1.size();
    // for (int i = 0; i< n; ++i) {
    //     if (NumericVector::is_na(nll1[i]))  nll1[i]=0;
    // }
    nll1[is_na(nll1)]=0;

    double nll2 = (0.79817986835+pow(fv,2))/2;

    return     sum(nll1) + nll2;
}
