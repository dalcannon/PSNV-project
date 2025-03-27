#include <Rcpp.h>
using namespace Rcpp;
// Function to double X
// [[Rcpp::export]]
NumericVector funcpp(NumericVector x) {
  return x * 2;
}

