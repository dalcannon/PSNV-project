#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix movmean(NumericMatrix X, int w) {
  int m=X.nrow(); // Number of rows in X
  int n=X.ncol(); // Number of columns in X
  Rcpp::NumericMatrix rollmean(m, n);
  if (w <= 0 || w % 2 == 0) {
    stop("w must be an odd positive integer.");
  }
  int v = (w - 1) / 2; // half-window size v = 2*w+1
  // Loop over each row in the matrix
  for (int i = 0; i < m; ++i) {
    // Loop over each column in the current row
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      int count = 0;
      for (int k = std::max(0, j - v); k <= std::min(n - 1, j + v); ++k) {
        sum += X(i, k);
        ++count;
      }
      rollmean(i,j) = sum/count;
    }
  }
  return wrap(rollmean);
}
