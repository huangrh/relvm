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
