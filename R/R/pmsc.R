#' Piecewise Multiplicative Scatter Correction (PMSC)
#'
#' @description
#' Perform piecewise multiplicative scatter correction (PMSC), i.e., MSC performed separately on a sliding band centered at the jth column of X
#'
#' @param X row vector or matrix of spectra to correct
#' @param w moving window size (must be an odd integer)
#' @param r reference spectrum
#' @param method Algorithm to calculate PSMC correction (see below for details)
#'
#' @return Matrix with PMSC-corrected spectra
#'
#' @details
#' 3 different algorithms for PMSC are supported using the 'methods' argument: "cpp" (default), "simple", and "mpmsc":
#'
#' "simple" will calculate the full MSC correction at a given moving window, and will extract the correction corresponding to the j-th wavelength.
#'
#' "cpp" calculates the MSC coefficients corresponding to the moving window surrounding the j-th column of X, and performs the correction using this manner.
#' This is all calculated inside a C++ wrapper function for speedup. This version is the fastest to calculate.
#'
#' "mpmsc" uses moving sums and moving means to calculate the correction; This version is optimal in MATLAB but is slightly slower in R
#'
#' @export
#' @examples
#' data("marzipan")
#' ref <- colMeans(marzipan$X)
#' X.pmsc <- mpmsc(marzipan$X,201,ref) #PMSC with window size of 201 channels
#' matplot(as.vector(marzipan$ax),t(X.pmsc),type="l",lty=1)
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) (Submitted)
#'
#' T. Isaksson, B. Kowalski, Applied Spectroscopy 1993, 47 (6), 702–709
#'
pmsc <- function( X, w, r, method=c("cpp","simple","mpmsc") ){
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
  if(missing(method)){method <- "cpp"}
  if(method %in% c("cpp","simple","mpmsc")==F){stop("Argument 'method' should be one of “cpp”, “simple”, “mpmsc”")}
  r <- matrix(r,nrow=1) #Ensure r is a ROW vector
  #Pre-compute intermediate moving results:
  switch(match.arg(method), simple = {
    X.pmsc <- X
    for(j in 1:n){
      band <- max(1,j-v):min(j+v,n)
      Xmsc <- MSC(X[,band],r = r[band])
      X.pmsc[,j] <- Xmsc[,which(band==j)]
    }
  }, cpp = {
    X.pmsc <- pmsc_cpp(X,w,r)
  }, mpmsc = {
    X.pmsc <- mpmsc(X,w,r)
  })
  return(X.pmsc)
}
