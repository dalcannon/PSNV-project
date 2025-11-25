#' Multiplicative Scatter Correction (MSC)
#'
#' @description
#' Performs multiplicative scatter correction (MSC). The code is a copy of the function 'msc' from the package 'pls', but has been slightly optimized for PMSC and LMSC
#' For standalone use, the code has no input sanitization or checks, so please ensure that class(X)="matrix" and that length(ref)=ncol(X) before use
#'
#' @param X matrix of spectra to correct (if X is a vector, use matrix(X,nrow=1) to convert to a row matrix)
#' @param r reference spectrum (default: average of each column of X)
#'
#' @return Matrix of MSC-corrected spectra
#'
#' @export
#' @examples
#' data("marzipan")
#' ref <- colMeans(marzipan$X)
#' X.msc <- MSC(marzipan$X,ref) #PMSC with window size of 101 channels
#' matplot(marzipan$ax,t(X.msc),type="l",lty=1)
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) (Submitted)
#' P Geladi, D MacDougall, H Martens, Applied Spectroscopy 1985, 39 (3), 491–500.

MSC <- function (X, r = colMeans(X)){
  Z <- cbind(1, r)
  B <- t(solve(crossprod(Z), t(X %*% Z)))
  res <- (X - B[, 1])/B[, 2]
  return(res)
}

