#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// Function to calculate the rolling standard deviation
// [[Rcpp::export]]
NumericMatrix movstd(Rcpp::NumericMatrix X, int w) {
  int m=X.nrow();
  int n=X.ncol();
  //Calculate the rolling mean (rollmean)
  if (w <= 0 || w % 2 == 0) {
    stop("w must be an odd positive integer.");
  }
  int v = (w - 1) / 2; // half-window size v = 2*w+1
  Rcpp::NumericMatrix rollmean(m, n);
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
  //now calculate the rolling standard deviation
  Rcpp::NumericMatrix movsd(m, n);
  //Rcpp::NumericMatrix rollmean = movmean(X, w); //calculate the rolling average of X
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double sumSquaredDiff = 0.0;
      int count = 0;
      for (int k = j - v; k <= j + v; ++k) {
        if (k >= 0 && k < n) {
          double diff = X(i,k) - rollmean(i,j);
          sumSquaredDiff += diff * diff;
          count++;
        }
      }
      double variance = count > 0 ? sumSquaredDiff / (count - 1) : 0.0;
      movsd(i,j) = count > 0 ? std::sqrt(variance) : NA_REAL;
    }
  }
  return movsd;
}
