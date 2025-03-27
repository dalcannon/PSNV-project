#' Piecewise Standard Normal Variate
#'
#' @description
#' Calculate piecewise standard normal variate correction (PSNV), i.e. SNV evaluated using a rolling average and rolling standard deviation.
#' For a window size w (an odd integer between 3 and 2:ncol(X)), the half-window size v is (w-1)/2. The moving window for channel j is equal to
#' max(1,j-v) to min(j+v,ncol(X)). As a result, at the edges of the window are asymmetric at the edges of the spectrum
#'
#' @param X vector or matrix of spectra to correct
#' @param w moving window size (must be an odd integer)
#'
#' @return Matrix with PSNV-corrected spectra
#'
#' @author Cannon Giglio
#' @export
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) [In Progress]
#' RJ Barnes, Mewa Singh Dhanoa, Susan J Lister, Applied Spectroscopy 1989, 43 (5), 772–777.
#' @examples
#' data(marzipan)
#' X.psnv <- psnv(marzipan$X,101) #PSNV with window size of 201 channels
#' matplot(as.vector(marzipan$ax),t(X.psnv),type="l",lty=1)
#'
#'
psnv <- function(X,w){
  if(!is.matrix(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  if(w%%2 != 1){stop("w must be an odd positive integer")}
  if(w<=1){stop("w too small")}
  if(w>(2*ncol(X))){stop("w too large")}
  v <- (w-1)/2
  rollparams <- movmeanandstd(X,v) #movmeanstd.cpp
  Xpsnv <- (X-rollparams$mean)/rollparams$sd
  return(Xpsnv)
}
