#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix pmsc_cpp(const NumericMatrix& Xr, int w, const NumericVector& rr) {
  const Map<MatrixXd> X(as<Map<MatrixXd>>(Xr)); // (m x n)
  const Map<VectorXd> r(as<Map<VectorXd>>(rr)); // length n

  int m = X.rows();
  int n = X.cols();
  if (w <= 0 || w % 2 == 0) {
    stop("w must be an odd positive integer.");
  }
  int v = (w - 1) / 2; // half-window size v = 2*w+1
  MatrixXd Xcorr = X; // output corrected data

  for (int j = 0; j < n; ++j) {
    int start = std::max(0, j - v);
    int end   = std::min(n - 1, j + v);
    int band_size = end - start + 1;

    // Subset X and r for the current window
    MatrixXd Xband = X.middleCols(start, band_size);
    VectorXd rband = r.segment(start, band_size);

    // Build Z = [1, rband]
    MatrixXd Z(band_size, 2);
    Z.col(0).setOnes();
    Z.col(1) = rband;

    // Compute pseudoinverse of Z and regression coefficients
    MatrixXd Zpinv = (Z.transpose() * Z).inverse() * Z.transpose(); // 2 x band_size
    MatrixXd B = Xband * Zpinv.transpose(); // m x 2

    // Correct column j: (X[,j] - B[,0]) / B[,1]
    for (int i = 0; i < m; ++i) {
      Xcorr(i, j) = (X(i, j) - B(i, 0)) / B(i, 1);
    }
  }

  return wrap(Xcorr);
}
