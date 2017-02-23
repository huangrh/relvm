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

