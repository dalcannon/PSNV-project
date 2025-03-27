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
#' X.pmsc <- pmsc(marzipan$X,101,ref) #PMSC with window size of 101 channels
#' matplot(marzipan$ax,t(X.pmsc),type="l",lty=1)
#'
pmsc2 <- function( X, w, r ){
  require(pls)
  #---------------------------------------------------------------------
  # PURPOSE: Perform piecewise MSC, i.e., MSC performed separately on
  # a sliding band centered at the jth wave(length/number)
  #---------------------------------------------------------------------
  #INPUTS:
  # [1] X: (m,n) matrix of spectra---each spectrum is aligned row-wise
  # [2] w: width of sliding band of wavelengths/wavenumbers (must be odd integer)
  # [3] r: n-vector---the reference spectrum
  #---------------------------------------------------------------------
  #OUTPUT:
  # [1] Xt: (m,n) matrix of transformed spectra
  #---------------------------------------------------------------------
  if(is.data.frame(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  if(missing(r)){r <- colMeans(X)}
  n <- ncol(X)# = size(X,2);
  if(length(r)!=n){stop("ref must match number of columns of X")}
  if(w%%2 != 1){stop("w must be an odd positive integer")}
  if(w>(2*n)){stop("w too large")}
  v <- (w-1)/2
  Xt <- X
  for(j in 1:n){
    band <- max(1,j-v):min(j+v,n)
    Xmsc <- pls::msc(X[,band],r[band])
    Xt[,j] <- Xmsc[,which(band==j)]
  }
}
