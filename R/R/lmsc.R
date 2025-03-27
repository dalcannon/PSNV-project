#' Localized Multiplicative Scatter Correction (LMSC)
#'
#' @description
#' Perform localized multiplicative scatter correction (LMSC), i.e., MSC performed separately on k disjoint but continuous spectral segments.
#'
#' @param X row vector or matrix of spectra to correct
#' @param k number of segments
#' @param r reference spectrum
#'
#' @return Matrix of LMSC-corrected spectra
#'
#' @export
#' @examples
#' data(marzipan)
#' ref <- colMeans(marzipan$X)
#' X.lmsc <- lmsc(marzipan$X,5,ref) #LMSC with 5 segments
#' matplot(as.vector(marzipan$ax),t(X.lmsc),type="l",lty=1)
#' @references C. Giglio, J.-M. Roger, E. Andries. Piecewise Standard Normal Variate (PSNV) for scatter correction. J Chemometrics (2025) [In Progress]
#' C.A. Andersson, Chemometrics and Intelligent Laboratory Systems 1999, 47 (1), 51–63.
#'
lmsc <- function(X,k,r){
  if(!is.matrix(X)){ X <- as.matrix(X)}
  if(is.null(nrow(X))||ncol(X)==1){ #x is a vector or a column matrix (1 x n)
    X <- matrix(X,nrow=1) #convert to a row matrix
  }
  if(missing(r)){r <- colMeans(X)}
  n <- ncol(X)
  if(length(r)!=n){stop("ref must match number of columns of X")}
  sidx <- c(0,which(diff(sort(((0:(n-1))%%k)+1))==1),n) #Segment indices
  Xnew <- X
  for(i in 1:k){
    band <- (sidx[i]+1):sidx[i+1]
    Xnew[,band] <- MSC(X[,band],r[band]) #fastmsc
  }
  return(Xnew)
}
