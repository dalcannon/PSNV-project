#' Piecewise Multiplicative Scatter Correction
#'
#' @description
#' Perform piecewise multiplicative scatter correction (PMSC), i.e., MSC performed separately on a sliding band centered at the jth wave(length/number)
#'
#' @param X row vector or matrix of spectra to correct
#' @param w moving window size (must be an odd integer)
#' @param r reference spectrum
#'
#' @return Matrix with PMSC-corrected spectra
#'
#' @export
#' @examples
#' data(marzipan)
#' ref <- colMeans(marzipan$X)
#' X.pmsc <- pmsc(marzipan$X,201,ref) #PMSC with window size of 201 channels
#' matplot(as.vector(marzipan$ax),t(X.pmsc),type="l",lty=1)
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) [In Progress]
#' T. Isaksson, B. Kowalski, Applied Spectroscopy 1993, 47 (6), 702–709
#'
pmsc <- function( X, w, r ){
  if(!is.matrix(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  if(missing(r)){r <- colMeans(X)}
  n <- ncol(X)# = size(X,2);
  if(length(r)!=n){stop("ref must match number of columns of X")}
  if(w%%2 != 1){stop("w must be an odd positive integer")}
  if(w<=1){stop("w too small")}
  if(w>(2*ncol(X))){stop("w too large")}
  v <- (w-1)/2
  Xpmsc <- X
  for(j in 1:n){
    band <- max(1,j-v):min(j+v,n)
    Xmsc <- MSC(X[,band],r[band])
    Xpmsc[,j] <- Xmsc[,which(band==j)]
  }
  return(Xpmsc)
}
